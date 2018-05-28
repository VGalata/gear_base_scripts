#!/bin/python
# -*- coding: utf-8 -*-

try:
    import _preamble
except ImportError:
    sys.exc_clear()

import os
from os import path
import sys
from sys import stdout
import pickle
import shutil
import logging
import re
import time, datetime

import itertools
from multiprocessing import Pool
import shlex, subprocess

from cv_class import CVClass
from utils.my_utils import CustomParser, timestamp_, timestamp, print_list
from utils.my_utils import make_str_compliant, run_cmd
from utils.cv_utils import arg_parser, add_args, get_p_f_s
from utils.cv_utils import cv_create_odirs, cv_samples_intersect, cv_pheno_stats, cv_pheno_check, cv_create_folds, cv_check_samples, cv_preproc_X_bin, cv_convert_X_vcf, cv_preproc_X_gds, cv_preproc_X_bin_vcf, cv_merge_Xs
from utils.cv_utils import cv_eig_convert, cv_eig_pca, cv_eig_post_proc, cv_eig_sel_k, cv_eig_full_post_proc
from utils.cv_utils import cv_model, cv_fold_sum, cv_rep_sum, cv_combi_sel_features, cv_clean, cv_pheno_sum


def write_log(info, logging, verbose=True):
    logging.info(info)
    if verbose: stdout.write(info); stdout.flush()

def save_cv_obj(cv_obj, cv_file, protocol=3, verbose=True):
    with open(cv_file, 'wb') as output:
        pickle.dump(cv_obj, output, protocol)

def load_cv_obj(cv_file, verbose=True):
    with open(cv_file, 'rb') as input:
        return( pickle.load(input) )

# NOTE MAIN
if __name__ == "__main__":
    #--------------------------------------------------#
    # Args
    parser = arg_parser()
    args = parser.parse_args()
    # Add some args
    args = add_args(args)
    print(args)
    #--------------------------------------------------#
    # Loggging
    logging.basicConfig(filename=args.log, level=logging.DEBUG)
    info ="\nBinary classification using k-fold CV on binary features\nStart: %s\n\n%s\n" % (timestamp(), parser.print_args(args))
    write_log(info, logging, args.verbose)
    
    #--------------------------------------------------#
    # X type is VCF -> GDS -> bin. matrix
    if args.X_type == 'vcf':
        args.X_gds = path.join(args.odir,'X.gds')
        args.X_bin = path.join(args.odir,'X.bin')
        vcf_gds_skip = ( path.isfile(args.X_gds) and path.isfile(args.X_bin) )
        info = 'Converting X VCF into GDS and BIN X\n'; write_log(info, logging, args.verbose)
        res = cv_convert_X_vcf(cmd=args.xvcf_convert_cmd, src_path=args.src_path, vcf_file=args.X_file, gds_file=args.X_gds, bin_file=args.X_bin, skip=vcf_gds_skip)
        assert res[2]==0
        info = "CMD: %s : %s\n%s\n" % (res[0],res[2],res[1]); write_log(info, logging, args.verbose)
    elif args.X_type == 'bin':
        args.X_bin = args.X_file
    elif args.X_type == 'bin_vcf':
        # some assertions
        assert args.X_file2 is not None
        assert len(args.proc_filter)==2
        assert len(args.proc_cmd)==2
        # VCF -> GDS -> Bin
        args.X_gds      = path.join(args.odir,'X.gds')
        args.X_bin      = path.join(args.odir,'X.bin')
        vcf_gds_skip = ( path.isfile(args.X_gds) and path.isfile(args.X_bin) )
        info = 'Converting X VCF into GDS and BIN X\n'; write_log(info, logging, args.verbose)
        res = cv_convert_X_vcf(cmd=args.xvcf_convert_cmd, src_path=args.src_path, vcf_file=args.X_file2, gds_file=args.X_gds, bin_file=args.X_bin, skip=vcf_gds_skip)
        assert res[2]==0
        info = "CMD: %s : %s\n%s\n" % (res[0],res[2],res[1]); write_log(info, logging, args.verbose)
        # VCF/GDS Bin + Bin
        info = 'Merging X VCF/GDS BIN and X BIN\n'; write_log(info, logging, args.verbose)
        # VCF/GDS BIN copy
        res = run_cmd("cp %s %s" % (args.X_bin,args.X_bin+'.tmp')) # make copy of VCF/GDS BIN but name ".tmp"
        assert res[1]==0
        merge_bin_skip = path.isfile(args.X_bin+'.gds') # skip if copy exists with ".gds"
        res = cv_merge_Xs(args.src_path, args.X_bin, args.in_samples, args.ex_samples, merge_bin_skip, args.X_file, args.X_bin)
        info = "CMD: %s : %s\n%s\n" % (res[0],res[2],res[1]); write_log(info, logging, args.verbose)
        assert res[2]==0
        res = run_cmd("mv %s %s" % (args.X_bin+'.tmp',args.X_bin+'.gds')) # rename copy
        assert res[1]==0
    else:
        sys.exit("Unknown X matrix type %s." % args.X_type)
    
    #--------------------------------------------------#
    features_f          = path.join(args.odir,"features.txt") # all features from X
    samples_f           = path.join(args.odir,"samples.txt")  # all samples from X
    
    # Phenotypes from Y, features and samples from X
    phenos, features, samples = get_p_f_s(args.Y_file, args.X_bin, x_s_in_header=True)
    info = "\nPhenos: [%d]\n%s\n" % (len(phenos),print_list(phenos,10)); write_log(info, logging, args.verbose)
    info = "\nFeatures in X: [%d]\n" % len(features); write_log(info, logging, args.verbose)
    info = "\nSamples in X: [%d]\n" % len(samples); write_log(info, logging, args.verbose)
    
    if not path.isfile(features_f):
        with open(features_f,'w') as f_f:
            for feature in features: f_f.write('%s\n' % feature)
    if not path.isfile(samples_f):
        with open(samples_f,'w') as f_f:
            for f_sample in samples: f_f.write('%s\n' % f_sample)
    
    #--------------------------------------------------#
    # Initiliaze CV object for each phenotype
    info = "\n[%s] Initializing CV\n" %(timestamp()); write_log(info, logging, args.verbose)
    if not path.exists("%s/obj" % args.odir):
        os.makedirs("%s/obj" % args.odir)
    
    cv_phenos = {}
    cv_obj    = {}
    x_files   = args.X_file if args.X_file2 is None else ','.join([args.X_file,args.X_file2]) # NEW
    for pheno in phenos:
        pheno_pickle = path.join("%s/obj" % args.odir,"%s.pkl" % make_str_compliant(pheno))
        if path.isfile(pheno_pickle):
            info = "Pheno %s: CV obj. exists: %s -> LOAD\n" % (pheno, pheno_pickle); write_log(info, logging, args.verbose)
            cv_phenos[pheno] = load_cv_obj(pheno_pickle)
            cv_obj[pheno] = pheno_pickle
        else:
            info = "Pheno %s: CV obj. does not exist: %s -> CREATE and SAVE\n" % (pheno, pheno_pickle); write_log(info, logging, args.verbose)
            cv_phenos[pheno] = CVClass( \
                X_file=x_files, X_type=args.X_type, X_bin=args.X_bin, X_gds=args.X_gds, X_vcf_add=None, \
                Y_file=args.Y_file, Y_type=args.Y_type, Y_name=pheno, \
                in_samples=args.in_samples, ex_samples=args.ex_samples, features = features_f, \
                odir=args.odir, \
                cv_reps=args.reps, cv_folds=args.folds \
            )
            cv_obj[pheno] = pheno_pickle
            save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        # add dicts
        cv_phenos[pheno].add_eig()
        cv_phenos[pheno].add_mod(args.model)
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
    
    #--------------------------------------------------#
    # Force steps:
    info=""
    for pheno, pheno_cv in cv_phenos.items():
        for f_s in args.force:
            if f_s in pheno_cv.done and (args.force_y is None or pheno in args.force_y):
                pheno_cv.done[f_s] = False #pheno_cv.done[f_s] #and not args.force[f_s]
                info += "FORCE: %s for %s\n" % (f_s,pheno)
            if f_s in pheno_cv.done[args.model] and (args.force_y is None or pheno in args.force_y):
                pheno_cv.done[args.model][f_s] = False #pheno_cv.done[args.model][f_s] #and not args.force[f_s]
                info += "FORCE: %s for %s\n" % (f_s,pheno)
    write_log(info+'\n', logging, args.verbose)
    
    #--------------------------------------------------#
    # CV next init.
    for pheno, pheno_cv in cv_phenos.items():
        # Create output dir.
        cv_create_odirs(pheno_cv)
    
    # Samples
    pool = Pool(args.cores)
    pool_iter = itertools.product(cv_phenos.values(),[args.src_path])
    results = pool.starmap( cv_samples_intersect , pool_iter )
    pool.close(); pool.join()
    for res in results:
        info = "\nSamples for %s:\n%s : %s\n%s\n" % (pheno_cv.Y_name,res[0],res[2],res[1]); write_log(info, logging, False)
        assert res[2]==0,"CMD status %s : %s" % (res[0],res[2])
    pheno_cv.done['samples'] = True
    save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
    
    for pheno, pheno_cv in cv_phenos.items():
        info = "\n%s\n[%s] Phenotype %s next init.\n" %("#%s#" % ("-"*50),timestamp(),pheno)
        write_log(info, logging, args.verbose)
        # Pheno class stats and check
        #if not pheno_cv.done['y_stats']:
        info = cv_pheno_stats(pheno_cv)
        write_log(info+'\n', logging, args.verbose)
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Pheno check
        #if not pheno_cv.done['y_check']:
        info = cv_pheno_check(pheno_cv,args.min_samples,args.max_class_ratio)
        write_log(info, logging, args.verbose)
        info = "\nPhenotype check for %s: %s\n" % (pheno_cv.Y_name,str(pheno_cv.Y_check)); write_log(info, logging, args.verbose)
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Print CV object
        info = "\n"+pheno_cv.print_cv_class(); write_log(info, logging, args.verbose)
        
        # Skip if not passed pheno check
        if not pheno_cv.Y_check:
            info = "\nSince check not passed skip (dir. removed)\n"
            write_log(info, logging, args.verbose)
            if path.isdir(pheno_cv.odir): shutil.rmtree(pheno_cv.odir)
            continue
        
        # Fold samples
        s = cv_create_folds(pheno_cv)
        info = "\nCV folds:\n%s\n" % s; write_log(info, logging, args.verbose)
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Sample check
        check_info = cv_check_samples(in_samples_f=args.in_samples, ex_samples_f=args.ex_samples, cv_obj=pheno_cv)
        info = "\nSample check for %s:\n%s\n" % (pheno_cv.Y_name,check_info); write_log(info, logging, args.verbose)
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
    
    #--------------------------------------------------#
    
    # For each pheno do work
    info = "\n\n\n[%s] Starting CV" %(timestamp())
    write_log(info, logging, args.verbose)
    #for pheno, pheno_cv in cv_phenos.items():
    for pheno in sorted(cv_phenos.keys()):
        pheno_cv = cv_phenos[pheno]
        if not pheno_cv.Y_check: continue
        
        info = "\n%s\n[%s] Phenotype %s START\n" %("#%s#" % ("-"*50),timestamp(),pheno)
        write_log(info, logging, args.verbose)
        
        # NOTE FULL GWAS
        # Pre-proc
        info = "\n[%s] FULL X PRE-PROC\n" %(timestamp()); write_log(info, logging, args.verbose)
        if args.X_type=='bin':
            res  = cv_preproc_X_bin(cv_obj=pheno_cv, cmd=args.proc_cmd[0], cmd_filter=args.proc_filter[0], src_path=args.src_path, rep=None, fold=None, cores=min(args.cores,5), skip=pheno_cv.done['full_preproc'])
            info = "FULL PRE-PROC:\n%s : %s\n%s\n" % (res[2],res[4],res[3]); write_log(info, logging, False)
            assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        elif args.X_type=='vcf':
            res  = cv_preproc_X_gds(cv_obj=pheno_cv, cmd=args.proc_cmd[0], cmd_filter=args.proc_filter[0], src_path=args.src_path, rep=None, fold=None, skip=pheno_cv.done['full_preproc'])
            info = "FULL PRE-PROC:\n%s : %s\n%s\n" % (res[2],res[4],res[3]); write_log(info, logging, False)
            assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        elif args.X_type=='bin_vcf': # TEST
            res = cv_preproc_X_bin_vcf(cv_obj=pheno_cv, cmd_bin=args.proc_cmd[0], cmd_filter_bin=args.proc_filter[0], cmd_gds=args.proc_cmd[1], cmd_filter_gds=args.proc_filter[1], src_path=args.src_path, rep=None, fold=None, cores=min(args.cores,5), skip=pheno_cv.done['full_preproc'])
            info = "FULL PRE-PROC:\n%s : %s\n%s\n" % (res[2],res[4],res[3]); write_log(info, logging, False)
            assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        pheno_cv.done['full_preproc'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        # Convert
        info = "\n[%s] FULL EIGENSTRAT CONVERT\n" %(timestamp()); write_log(info, logging, args.verbose)
        skip_convert = pheno_cv.done['full_eig_conv'] or (pheno_cv.done['full_eig_pca'] and pheno_cv.done['full_eig_assoc'] and pheno_cv.done['full_eig_proc'])
        cv_eig_convert(cv_obj=pheno_cv, rep=None, fold=None, skip=skip_convert)
        pheno_cv.done['full_eig_conv'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        # PCA
        info = "\n[%s] FULL EIGENSTRAT PCA\n" %(timestamp()); write_log(info, logging, args.verbose)
        res  = cv_eig_pca(cv_obj=pheno_cv, eig_path=args.eig_path, cmd=args.eig_pca_cmd, k=args.eig_k_max, rep=None, fold=None, skip=pheno_cv.done['full_eig_pca'])
        info = "FULL EIGENSTRAT PCA:\n%s : %s\n%s\n" % (res[2],res[4],res[3]); write_log(info, logging, False)
        assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        pheno_cv.done['full_eig_pca'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        # Assoc + select k w.r.t. lambda
        info = "\n[%s] FULL EIGENSTRAT ASSOC\n" %(timestamp()); write_log(info, logging, args.verbose)
        res  = cv_eig_sel_k(pheno_cv,args.eig_path,args.eig_assoc_cmd,args.eig_lambda_cmd,k_max=args.eig_k_max,k_step=args.eig_k_step,l_min=args.eig_l_min, rep=None, fold=None, skip=pheno_cv.done['full_eig_assoc'])
        info = "FULL EIGENSTRAT ASSOC:\n%s\n" % (res[2]); write_log(info, logging, False)
        pheno_cv.done['full_eig_assoc'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        # Post-processing
        info = "\n[%s] FULL EIGENSTRAT POST-PROC\n" %(timestamp()); write_log(info, logging, args.verbose)
        res  = cv_eig_full_post_proc(cv_obj=pheno_cv, src_path=args.src_path, cmd=args.eig_proc_full_cmd, cores=min(args.cores,30), skip=pheno_cv.done['full_eig_proc'])
        info = "FULL EIGENSTRAT POST-PROC:\n%s : %s\n%s\n" % (res[0],res[2],res[1]); write_log(info, logging, False)
        assert res[2]==0,"CMD status %s : %s" % (res[2],res[0])
        pheno_cv.done['full_eig_proc'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # NOTE REPS/FOLDS EIGENSTRAT
        # Pre-processing
        info = "\n[%s] X PRE-PROC\n" %(timestamp()); write_log(info, logging, args.verbose)
        pool = Pool(max(1,int(args.cores/3))) # NOTE changed from 2 to 3 to have less concurrent jobs because R-package parallel's parSapply in my_utils.R may freeze
        if args.X_type=='bin':
            # cv_preproc_X_bin(cv_obj, cmd, cmd_filter, src_path, rep=None, fold=None, cores=2, skip=False)
            # NOTE you may want to set cores (second last param) to 1 if pipeline is getting stuck here every time
            pool_iter = itertools.product([pheno_cv],[args.proc_cmd[0]],[args.proc_filter[0]],[args.src_path],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[2],[pheno_cv.done['cv_preproc']])
            results = pool.starmap( cv_preproc_X_bin , pool_iter )
        elif args.X_type=='vcf':
            pool_iter = itertools.product([pheno_cv],[args.proc_cmd[0]],[args.proc_filter[0]],[args.src_path],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[pheno_cv.done['cv_preproc']])
            results = pool.starmap( cv_preproc_X_gds , pool_iter )
        elif args.X_type=='bin_vcf': # TEST
            pool_iter = itertools.product([pheno_cv],[args.proc_cmd[0]],[args.proc_filter[0]],[args.proc_cmd[1]],[args.proc_filter[1]],[args.src_path],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[2],[pheno_cv.done['cv_preproc']])
            results = pool.starmap( cv_preproc_X_bin_vcf , pool_iter )
        pool.close(); pool.join()
        for res in results:
            info="rep %d, fold%d, X PRE-PROC:\n%s : %s\n%s\n" % (res[0],res[1],res[2],res[4],res[3]); write_log(info, logging, False)
            assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        pheno_cv.done['cv_preproc'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Convert
        info = "\n[%s] EIGENSTRAT CONVERT\n" %(timestamp()); write_log(info, logging, args.verbose)
        skip_convert = pheno_cv.done['cv_eig_conv'] or (pheno_cv.done['cv_eig_pca'] and pheno_cv.done['cv_eig_assoc'] and pheno_cv.done['cv_eig_proc'])
        pool = Pool(args.cores)
        pool_iter = itertools.product([pheno_cv],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[skip_convert])
        results = pool.starmap( cv_eig_convert , pool_iter )
        pool.close(); pool.join()
        pheno_cv.done['cv_eig_conv'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # PCA: cv_obj, eig_path, cmd, k=10, rep=None, fold=None, skip=False
        info = "\n[%s] EIGENSTRAT PCA\n" %(timestamp()); write_log(info, logging, args.verbose)
        pool = Pool(max(1,int(args.cores/2))) # NOTE smartpca uses multithreading but only used for a very short time so
        pool_iter = itertools.product([pheno_cv],[args.eig_path],[args.eig_pca_cmd],[args.eig_k_max],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[pheno_cv.done['cv_eig_pca']])
        results = pool.starmap( cv_eig_pca , pool_iter )
        pool.close(); pool.join()
        for res in results:
            info="rep %d, fold%d, EIGENSTRAT PCA:\n%s : %s\n%s\n" % (res[0],res[1],res[2],res[4],res[3]); write_log(info, logging, False)
            assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        pheno_cv.done['cv_eig_pca'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # ASSOC: cv_eig_sel_k(cv_obj, eig_path, eig_assoc_cmd, eig_lambda_cmd, k_max, l_max=1.1, rep=None, fold=None, skip=False)
        info = "\n[%s] EIGENSTRAT ASSOC\n" %(timestamp()); write_log(info, logging, args.verbose)
        pool = Pool(args.cores)
        pool_iter = itertools.product([pheno_cv],[args.eig_path],[args.eig_assoc_cmd],[args.eig_lambda_cmd],[args.eig_k_max],[args.eig_k_step],[args.eig_l_min],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[pheno_cv.done['cv_eig_assoc']])
        results = pool.starmap( cv_eig_sel_k , pool_iter )
        pool.close(); pool.join()
        for res in results:
            info="rep %d, fold%d, EIGENSTRAT ASSOC:\n%s\n" % (res[0],res[1],res[2])
            write_log(info, logging, False)
        pheno_cv.done['cv_eig_assoc'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Post-processing
        info = "\n[%s] EIGENSTRAT POST-PROC\n" %(timestamp()); write_log(info, logging, args.verbose)
        pool = Pool(args.cores)
        pool_iter = itertools.product([pheno_cv],[args.src_path],[args.eig_proc_cmd],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[pheno_cv.done['cv_eig_proc']])
        results = pool.starmap( cv_eig_post_proc , pool_iter )
        pool.close(); pool.join()
        for res in results:
            info="rep %d, fold%d, EIGENSTRAT POST-PROC:\n%s : %s\n%s\n" % (res[0],res[1],res[2],res[4],res[3]); write_log(info, logging, False)
            assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        pheno_cv.done['cv_eig_proc'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # MODELS
        info = "\n[%s] MODELS\n" %(timestamp()); write_log(info, logging, args.verbose)
        pool = Pool(args.cores)
        pool_iter = itertools.product([pheno_cv],[args.src_path],[args.model],[args.model_cmd],[args.model_params],range(1,pheno_cv.reps+1),range(1,pheno_cv.folds+1),[pheno_cv.done[args.model]['cv_mod']])
        results = pool.starmap( cv_model , pool_iter )
        pool.close(); pool.join()
        for res in results:
            info="rep %d, fold%d, MODELS:\n%s : %s\n%s\n" % (res[0],res[1],res[2],res[4],res[3]); write_log(info, logging, False)
            assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        pheno_cv.done[args.model]['cv_mod'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Fold summary
        info = "\n[%s] FOLD SUM\n" %(timestamp()); write_log(info, logging, args.verbose)
        pool = Pool(args.cores)
        pool_iter = itertools.product([pheno_cv],[args.src_path],[args.model],range(1,pheno_cv.reps+1),[pheno_cv.done[args.model]['cv_fold_sum']])
        results = pool.starmap( cv_fold_sum , pool_iter )
        pool.close(); pool.join()
        for res in results:
            info = "rep %d, FOLD SUM:\n%s : %s\n%s\n" % (res[0],res[1],res[3],res[2])
            write_log(info, logging, False)
            assert res[3]==0,"CMD status %s : %s" % (res[3],res[2])
        pheno_cv.done[args.model]['cv_fold_sum'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Rep summary
        info = "\n[%s] REP SUM\n" %(timestamp()); write_log(info, logging, args.verbose)
        res  = cv_rep_sum(cv_obj=pheno_cv, src_path=args.src_path, model=args.model, skip=pheno_cv.done[args.model]['cv_rep_sum'])
        info ="FOLD SUM:\n%s : %s\n%s\n" % (res[0],res[2],res[1]); write_log(info, logging, False)
        assert res[2]==0,"CMD status %s : %s" % (res[2],res[0])
        pheno_cv.done[args.model]['cv_rep_sum'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # NOTE FULL
        # Full: selected features
        info = "\n[%s] SEL FEATURES FOR FULL\n" %(timestamp()); write_log(info, logging, args.verbose)
        res  = cv_combi_sel_features(cv_obj=pheno_cv, src_path=args.src_path, skip=pheno_cv.done['full_sel_f'])
        info = "SEL FEATURES FOR FULL:\n%s : %s\n%s\n" % (res[0],res[2],res[1]); write_log(info, logging, False)
        assert res[2]==0,"CMD status %s : %s" % (res[2],res[0])
        pheno_cv.done['full_sel_f'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # Model
        info = "\n[%s] FULL MODEL\n" %(timestamp()); write_log(info, logging, args.verbose)
        res  = cv_model(pheno_cv, src_path=args.src_path, model=args.model, model_cmd=args.model_full_cmd, model_params=args.model_params, rep=None, fold=None, skip=pheno_cv.done[args.model]['full_mod'])
        info = "FULL MODEL:\n%s : %s\n%s\n" % (res[2],res[4],res[3]); write_log(info, logging, False)
        assert res[4]==0,"CMD status %s : %s" % (res[4],res[2])
        pheno_cv.done[args.model]['full_mod'] = True
        save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        # NOTE Clean
        if args.clean_up:
            info = "\n[%s] CLEAN\n" %(timestamp()); write_log(info, logging, args.verbose)
            res = cv_clean(pheno_cv)
            info = "\n%s" % res; write_log(info, logging, False)
            save_cv_obj(cv_phenos[pheno],cv_obj[pheno])
        
        info = "\n[%s] Phenotype %s END\n%s\n" %(timestamp(),pheno,"#%s#" % ("-"*50))
        write_log(info, logging, args.verbose)
    #--------------------------------------------------#
    
    # Summary over phenos
    # Models: cv_objs, src_path, odir, skip=False
    info = "\n[%s] PHENO SUM\n" %(timestamp()); write_log(info, logging, args.verbose)
    res  = cv_pheno_sum(cv_objs=cv_phenos.values(), src_path=args.src_path, model=args.model, odir=args.odir)
    info = "PHENO SUM:\n%s : %s\n%s\n" % (res[0],res[2],res[1]); write_log(info, logging, False)
    