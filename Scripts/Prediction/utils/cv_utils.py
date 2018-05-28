#!/bin/python
# -*- coding: utf-8 -*-

# NOTE imports
import os
from os import path
import sys
from sys import stdout
import time
import re
import argparse
from sklearn import cross_validation
import numpy
import itertools
from multiprocessing import Pool
import shlex, subprocess
import shutil

from utils.my_utils import run_cmd, CustomParser, timestamp_, read_id_list, pheno_as_list, pheno_class_stats
from utils.my_utils import get_phenos, get_features, get_samples
from utils.convert_bin import bin_to_eig_geno_snp, bin_to_eig_pheno


# NOTE global vars
#step_choices = [\
    #'init_vcf','init_init','step_init',\
    #'cv_folds','cv_preproc','cv_eig_conv','cv_eig_pca','cv_eig_assoc','cv_eig_proc','cv_mod','cv_fold_sum','cv_rep_sum','step_cv',\
    #'full_sel_f','full_preproc','full_eig_conv','full_eig_pca','full_eig_tw','full_eig_lambda','full_eig_assoc','full_eig_proc','full_mod','step_full',\
    #'step_final'\
#]
step_choices = [\
    'samples','sample_check','y_stats','y_check'\
    'cv_folds','cv_preproc','cv_eig_conv','cv_eig_pca','cv_eig_assoc','cv_eig_proc',\
    'cv_mod','cv_fold_sum','cv_rep_sum',\
    'full_preproc','full_eig_conv','full_eig_pca','full_eig_assoc','full_eig_proc',
    'full_sel_f','full_mod'\
]

rm_files_eig = [\
    re.compile('.*/eig\.geno$'),\
    #re.compile('.*/eig\.snps$'),\
    re.compile('.*/eig\.pheno$')\
]
rm_files = []
rm_files.extend(rm_files_eig)


# NOTE methods
def arg_parser():
    # NOTE takes a file with prefix "@" as input for arguments
    parser = CustomParser(fromfile_prefix_chars='@', formatter_class=argparse.RawTextHelpFormatter)
    
    # NOTE Required arguments: input/output files/dir.s
    parser.add_argument('-X_file',                      help='Feature file (VCF or bin. features x samples); bin. if bin_vcf', required=True)
    parser.add_argument('-X_file2',                     help='VCF Feature file if bin_vcf')
    parser.add_argument('-Y_file',                      help='Phenotype file (samples x features)', required=True)
    parser.add_argument('-X_type',                      help='', choices=['bin','vcf','bin_vcf'], required=True)
    parser.add_argument('-Y_type',                      help='', choices=['bin'], default='bin')
    parser.add_argument('-odir',                        help='Output directory', required=True)
    
    # NOTE CV arguments
    parser.add_argument('--reps',                       help='Number of CV repetitions,', type=int, default=1)
    parser.add_argument('--folds',                      help='Number of folds in CV', type=int, default=5)
    parser.add_argument('--min_samples',                help='Minimal number of samples to run CV (after sample processing)', default=50, type=int)
    parser.add_argument('--max_class_ratio',            help='Maximal class ratio to run CV', default=75.0, type=float)
    
    # NOTE CMDs
    parser.add_argument('-proc_filter',                 help='if bin_vcf, first bin than VCF', required=True, nargs='+')
    parser.add_argument('-proc_cmd',                    help='if bin_vcf, first bin than VCF', required=True, nargs='+')
    
    parser.add_argument('-eig_path',                    help='', required=True)
    parser.add_argument('-eig_pca_cmd',                 help='', required=True)
    parser.add_argument('-eig_k_max',                   help='', default=10,    type=int)
    parser.add_argument('-eig_k_step',                  help='', default=1,     type=int)
    parser.add_argument('-eig_l_min',                   help='', default=1.0,   type=float)
    parser.add_argument('-eig_alpha',                   help='', default=0.05,  type=float)
    parser.add_argument('-eig_tw_cmd',                  help='')
    parser.add_argument('-eig_assoc_cmd',               help='', required=True)
    parser.add_argument('-eig_lambda_cmd',              help='', required=True)
    parser.add_argument('-eig_proc_cmd',                help='', required=True)
    parser.add_argument('-eig_proc_full_cmd',           help='', required=True)
    
    parser.add_argument('-model',                       help='', required=True, choices=['ctree','rpart'])
    parser.add_argument('-model_params',                help='', required=True)
    parser.add_argument('-model_cmd',                   help='', required=True)
    parser.add_argument('-model_full_cmd',              help='', required=True)
    
    ## NOTE Preprocessing VCF/GDS
    parser.add_argument('--xvcf_convert_cmd',           help='')
    
    ## NOTE Samples
    parser.add_argument('--in_samples',                 help='Sample list (one sample per line')
    parser.add_argument('--ex_samples',                 help='Sample list (one sample per line')
    
    # NOTE Other
    parser.add_argument('--src_path',                   help='', default='.')
    parser.add_argument('--clean_up',                   help='Will remove some intermediate files', action='store_true')
    parser.add_argument('--skip',                       help='Which steps to skip', nargs='+', choices=step_choices)
    parser.add_argument('--force',                      help='Which steps to force', nargs='+', choices=step_choices)
    parser.add_argument('--skip_y',                     help='For which phenotypes to skip steps', nargs='+')
    parser.add_argument('--force_y',                    help='For which phenotypes to force steps', nargs='+')
    parser.add_argument('--cores',                      help='Number of cores to use,', type=int, default=1)
    parser.add_argument('--dry_run',                    help='Do not execute CMDs', action='store_true')
    parser.add_argument('--verbose',                    help='Print additional information', action='store_true')
    
    return(parser)

def add_args(args):
    vars(args)['log'] = path.join(args.odir,"cv_%s.log" % timestamp_())
    #vars(args)['log'] = path.join(args.odir,"cv.log")
    vars(args)['X_bin'] = None
    vars(args)['X_gds'] = None
    #vars(args)['X_vcf_add'] = None
    if args.force is None: args.force=[]
    if args.skip is None: args.skip=[]
    return(args)

def get_p_f_s(Y_file, X_file, x_s_in_header=True):
    phenos      = get_phenos(Y_file)
    features    = get_features(X_file)
    samples     = get_samples(X_file, in_header=x_s_in_header)
    return phenos, features, samples
    

# NOTE Methods for CV object
"""
Create CV output directories
"""
def cv_create_odirs(cv_obj):
    if not path.exists(cv_obj.odir):
        os.makedirs(cv_obj.odir)
    if not path.exists(cv_obj.full['odir']):
        os.makedirs(cv_obj.full['odir'])
    for rep in cv_obj.cv.keys():
        for fold in cv_obj.cv[rep].keys():
            if not path.exists(cv_obj.cv[rep][fold]['odir']):
                os.makedirs(cv_obj.cv[rep][fold]['odir'])

"""
Sample intersection w.r.t. given sample lists
"""
def cv_samples_intersect(cv_obj, src_path):
    # already done
    if cv_obj.done['samples']:
        return ("SKIP: Sample file already created","",0)
    cmd = "Rscript {src_path}/utils/process_x_y_samples.R -geno_file {g} -pheno_file {p} -pheno_name {n} -o_samples {o} -o_stats {s} --rm_miss --miss_value NA --src_path {src_path}/utils --include_full --verbose"
    #cmd = cmd.format(src_path=src_path, p=cv_obj.Y_file, n=cv_obj.Y_name, o=cv_obj.samples, s=cv_obj.Y_stat['file'])
    cmd = cmd.format(src_path=src_path, p=cv_obj.Y_file, g=cv_obj.X_bin, n=cv_obj.Y_name, o=cv_obj.samples, s=cv_obj.Y_stat['file'])
    if cv_obj.in_samples is not None: cmd += " --include_samples %s" % cv_obj.in_samples
    if cv_obj.ex_samples is not None: cmd += " --exclude_samples %s" % cv_obj.ex_samples
    p_stdout, p_status = run_cmd(cmd)
    return (cmd, p_stdout, p_status)

"""
Pheno stat.s
"""
def cv_pheno_stats(cv_obj):
    # already done
    if cv_obj.done['y_stats']:
        return "SKIP: Pheno stat.s already computed"
    cv_obj.Y_stat['class_num'] = { '0' : 0 , '1' : 0 }
    cv_obj.Y_stat['class_pct'] = { '0' : 0 , '1' : 0 }
    with open(cv_obj.Y_stat['file']) as infile:
        for line in infile:
            line = line.rstrip("\n"); line = line.split('\t')
            if line[0]=='miss':
                cv_obj.Y_stat['miss'] = int(line[1])
            elif line[0]=='0':
                cv_obj.Y_stat['class_num']['0'] = int(line[1])
            elif line[0]=='1':
                cv_obj.Y_stat['class_num']['1'] = int(line[1])
            else:
                sys.exit("In Y stat.s file: unknown value %s" % line[0])
    cv_obj.Y_stat['total'] = sum([cv_obj.Y_stat['miss'],cv_obj.Y_stat['class_num']['0'],cv_obj.Y_stat['class_num']['1']])
    samples = set(read_id_list(cv_obj.samples))
    assert len(samples) == cv_obj.Y_stat['total'], "Have %d samples in %s but only %d as total from Y stat.s file %s" %(len(samples),cv_obj.samples,cv_obj.Y_stat['total'],cv_obj.Y_stat['file'])
    assert cv_obj.Y_stat['miss']==0, "Missing pheno for %s" % cv_obj.Y_name
    if (cv_obj.Y_stat['total']-cv_obj.Y_stat['miss']) > 0:
        cv_obj.Y_stat['class_pct'] = { \
            '0' : 100.0 * float(cv_obj.Y_stat['class_num']['0']) / float(cv_obj.Y_stat['total']-cv_obj.Y_stat['miss']) ,\
            '1' : 100.0 * float(cv_obj.Y_stat['class_num']['1']) / float(cv_obj.Y_stat['total']-cv_obj.Y_stat['miss']) \
        }
    cv_obj.done['y_stats'] = True
    return "Computed pheno stat.s"

"""
Pheno stat.s check: whether prediction can/should be done
"""
def cv_pheno_check(cv_obj, min_samples, max_class_ratio):
    if cv_obj.done['y_check']:
        return "SKIP: Pheno check already done"
    assert cv_obj.Y_stat['total'] is not None
    assert cv_obj.Y_stat['miss']  is not None
    assert cv_obj.Y_stat['class_pct'] is not None
    cv_obj.Y_check = (cv_obj.Y_stat['total']-cv_obj.Y_stat['miss']) >= min_samples and max(cv_obj.Y_stat['class_pct'].values()) <= max_class_ratio
    cv_obj.done['y_check'] = True
    return "Pheno check done"

"""
Create folds
"""
def cv_create_folds(cv_obj):
    if cv_obj.done['cv_folds']:
        return "SKIP: Folds already created"
    info = ""
    # did not pass the phenotype check
    if not cv_obj.Y_check:
        return "CV folds: %s did not pass pheno. check" % cv_obj.Y_name
    # Samples to use
    samples = set(read_id_list(cv_obj.samples))
    # Y as array and sample IDs
    y, y_samples = pheno_as_list(pheno_file=cv_obj.Y_file, pheno_name=cv_obj.Y_name, ignore_miss=True, samples=samples)
    y            = numpy.array(y)
    y_samples    = numpy.array(y_samples) # need for indexing
    for rep in cv_obj.cv.keys():
        # Create folds
        y_folds = cross_validation.StratifiedKFold(y=y, n_folds=cv_obj.folds, shuffle=True, random_state=None)
        # Save folds as lists of samples in train/test
        fold = 1
        for train_index, test_index in y_folds:
            with open(cv_obj.cv[rep][fold]['samples_train'],'w') as o_file:
                o_file.write("\n".join( y_samples[train_index].tolist() )+"\n")
            with open(cv_obj.cv[rep][fold]['samples_test'],'w') as o_file:
                o_file.write("\n".join( y_samples[test_index].tolist() )+"\n")
            # Stats
            y_fold_stats = pheno_class_stats(pheno_file=cv_obj.Y_file, pheno_name=cv_obj.Y_name, samples=set(y_samples[train_index].tolist()))
            info +=  "Rep %d, fold %d, train: %s\n" % (rep,fold,' ; '.join([' - '.join([k,"%d" % v]) for k,v in y_fold_stats.items()]))
            y_fold_stats = pheno_class_stats(pheno_file=cv_obj.Y_file, pheno_name=cv_obj.Y_name, samples=set(y_samples[test_index].tolist()))
            info +=  "Rep %d, fold %d, test: %s\n"  % (rep,fold,' ; '.join([' - '.join([k,"%d" % v]) for k,v in y_fold_stats.items()]))
            fold += 1
    cv_obj.done['cv_folds'] = True
    return info

"""
(Pre)processing of bin. feature matrix X
"""
def cv_preproc_X_bin(cv_obj, cmd, cmd_filter, src_path, rep=None, fold=None, cores=2, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    assert cv_obj.X_type == "bin", "Expected X type \"bin\" but have %s" % cv_obj.X_type
    o_f = None; in_s = None
    if rep is None or fold is None:
        o_f  = cv_obj.full['features_pr']
        in_s = cv_obj.full['samples']
    else:
        o_f  = cv_obj.cv[rep][fold]['features_pr']
        in_s = cv_obj.cv[rep][fold]['samples_train']
    cmd = cmd.format(src_path=src_path, mat=cv_obj.X_bin, o_f=o_f, in_s=in_s, proc_filter=cmd_filter, cores=cores)
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (rep, fold, cmd, p_stdout, p_status)

def cv_convert_X_vcf(cmd, src_path, vcf_file, gds_file, bin_file, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    cmd = cmd.format(src_path=src_path, v_file=vcf_file, g_file=gds_file, b_file=bin_file)
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (cmd, p_stdout, p_status)

def cv_preproc_X_gds(cv_obj, cmd, cmd_filter, src_path, rep=None, fold=None, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    assert cv_obj.X_type == "vcf", "Expected X type \"vcf\" but have %s" % cv_obj.X_type
    o_f = None; in_s = None
    if rep is None or fold is None:
        o_f  = cv_obj.full['features_pr']
        in_s = cv_obj.full['samples']
    else:
        o_f  = cv_obj.cv[rep][fold]['features_pr']
        in_s = cv_obj.cv[rep][fold]['samples_train']
    cmd = cmd.format(src_path=src_path, gds_file=cv_obj.X_gds, o_f=o_f, in_s=in_s, proc_filter=cmd_filter)
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (rep, fold, cmd, p_stdout, p_status)

def cv_preproc_X_bin_vcf(cv_obj, cmd_bin, cmd_filter_bin, cmd_gds, cmd_filter_gds, src_path, rep=None, fold=None, cores=2, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0; cmd = ''
    assert cv_obj.X_type == "bin_vcf", "Expected X type \"bin_vcf\" but have %s" % cv_obj.X_type
    o_f = None; in_s = None
    if rep is None or fold is None:
        o_f  = cv_obj.full['features_pr']
        in_s = cv_obj.full['samples']
    else:
        o_f  = cv_obj.cv[rep][fold]['features_pr']
        in_s = cv_obj.cv[rep][fold]['samples_train']
    o_f_1 = o_f + ".vcf"
    o_f_2 = o_f + ".bin"
    # VCF/GDS preproc
    cmd_gds = cmd_gds.format(src_path=src_path, gds_file=cv_obj.X_gds, o_f=o_f_1, in_s=in_s, proc_filter=cmd_filter_gds)
    if not skip:
        p_so, p_status = run_cmd(cmd_gds)
        p_stdout += ("\n" + p_so)
        cmd = cmd_gds
    if p_status!=0: # non-zero status -> return: calling code should check the exit status
        return (rep, fold, cmd, p_stdout, p_status)
    # Bin. preproc
    cmd_bin = cmd_bin.format(src_path=src_path, mat=cv_obj.X_file.split(',')[0], o_f=o_f_2, in_s=in_s, proc_filter=cmd_filter_bin, cores=cores)
    if not skip:
        p_so, p_status = run_cmd(cmd_bin)
        p_stdout += ("\n" + p_so)
        cmd += ("\n" + cmd_bin)
    if p_status!=0: # non-zero status -> return: calling code should check the exit status
        return (rep, fold, cmd, p_stdout, p_status)
    # Concat both feature lists
    cmd_cat = "cat %s %s" % (o_f_1,o_f_2)
    if not skip:
        p_so, p_status = run_cmd(cmd_cat)
        with open(o_f, "w") as f:
            p = subprocess.Popen(shlex.split(cmd_cat), stdout=f, stderr=subprocess.PIPE)
        p_so = p.stderr.read().decode()
        p_comm   = p.communicate()[0]
        p_status = p.returncode
        p_stdout += ("\n" + p_so)
        cmd += ("\n" + cmd_cat)
    return (rep, fold, cmd, p_stdout, p_status)

"""
"""
def cv_merge_Xs(src_path, ofile, in_samples=None, ex_samples=None, skip=False, *mat_files):
    p_stdout = 'SKIPPED'; p_status = 0
    cmd = "Rscript {src_path}/utils/merge_tables.R -mat_files {mat_files} -ofile {ofile} --verbose --src_path {src_path}/utils"
    cmd = cmd.format(src_path=src_path, mat_files=' '.join(mat_files), ofile=ofile)
    if in_samples is not None:
        cmd += (" --include_samples %s" % in_samples)
    if ex_samples is not None:
        cmd += (" --exclude_samples %s" % ex_samples)
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (cmd, p_stdout, p_status)

"""
"""
def cv_eig_convert(cv_obj, rep=None, fold=None, skip=False):
    eig_geno = eig_snps = eig_pheno = eig_snps_rm = None
    samples = features = features_pr = None
    cv_obj_sub = None
    if rep is None or fold is None:
        cv_obj_sub = cv_obj.full
    else:
        cv_obj_sub = cv_obj.cv[rep][fold]
    eig_geno     = cv_obj_sub['eig']['eig.geno']
    eig_snps     = cv_obj_sub['eig']['eig.snps']
    eig_snps_rm  = cv_obj_sub['eig']['eig.snps.rm']
    eig_pheno    = cv_obj_sub['eig']['eig.pheno']
    if rep is None or fold is None:
        samples  = cv_obj_sub['samples']
    else:
        samples  = cv_obj_sub['samples_train']
    features     = cv_obj_sub['features']
    features_pr  = cv_obj_sub['features_pr']
    if not skip:
        bin_to_eig_geno_snp( \
            geno_file=cv_obj.X_bin, \
            geno_ofile=eig_geno, snp_ofile=eig_snps, \
            sample_file=samples, feature_file=features, \
            miss='NA' \
        )
        bin_to_eig_pheno( \
            pheno_file=cv_obj.Y_file, \
            pheno_name=cv_obj.Y_name, \
            pheno_ofile=eig_pheno, \
            sample_file=samples, \
            miss='NA' \
        )
        snps    = set(read_id_list(features))
        snps_pr = set(read_id_list(features_pr))
        snps_rm = snps.difference(snps_pr)
        with open(eig_snps_rm,'w') as ofile:
            for snp in snps_rm: ofile.write("rs%s\n" % snp)

"""
"""
def cv_eig_pca(cv_obj, eig_path, cmd, k=10, rep=None, fold=None, skip=False):
    eig_pca = eig_ev = eig_pca_log = eig_pca_plot = None
    p_stdout = 'SKIPPED'; p_status = 0
    cv_obj_sub = None
    if rep is None or fold is None:
        cv_obj_sub = cv_obj.full
    else:
        cv_obj_sub = cv_obj.cv[rep][fold]
    assert 'eig' in cv_obj_sub
    eig_pca     = cv_obj_sub['eig']['eig.pca']
    eig_ev      = cv_obj_sub['eig']['eig.ev']
    eig_pca_log = cv_obj_sub['eig']['eig.pca.log']
    eig_pca_plot= cv_obj_sub['eig']['eig.pca.plot']
    eig_geno    = cv_obj_sub['eig']['eig.geno']
    eig_snps    = cv_obj_sub['eig']['eig.snps']
    eig_snps_rm = cv_obj_sub['eig']['eig.snps.rm']
    eig_pheno   = cv_obj_sub['eig']['eig.pheno']
    
    cmd = cmd.format( \
        eig_path=eig_path, eig_k=k, \
        eig_geno=eig_geno, eig_snps=eig_snps, eig_snps_rm=eig_snps_rm, eig_pheno=eig_pheno, \
        eig_pca=eig_pca,   eig_ev=eig_ev,     eig_pca_log=eig_pca_log, eig_pca_plot=eig_pca_plot \
    )
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (rep, fold, cmd, p_stdout, p_status)

"""
{eig_path}/bin/twstats -t {eig_path}/POPGEN/twtable -i {eig_ev} -o {eig_pca_pv}
"""
def cv_eig_twstats_k(eig_pca_pv, k_max=10, alpha=0.05):
    k = pv_co = n = 0
    header = True
    with open(eig_pca_pv,'r') as pc_pv:
        for line in pc_pv:
            line = line.rstrip("\n")
            line = line.split(' ')
            line = [l for l in line if l!=''] # remove empty
            if header:
                header = False
                n = len(line) - 1 # "effect .n" is split into two fields
                for i in range(0,len(line)):
                    if line[i] == 'p-value':
                        pv_col = i
                        break
            assert len(line) == n, "In %s: %s" % (eig_pca_pv,' '.join(line))
            if float(line[pv_col]) <= alpha : k += 1
            else: break
    return min([max([1,k]),k_max])

"""
"""
def cv_eig_twstats(cv_obj, eig_path, cmd, rep=None, fold=None, k_max=10, alpha=0.05, skip=False):
    eig_pca_pv = None; eig_k = None
    p_stdout = 'SKIPPED'; p_status = 0
    cv_obj_sub = None
    if rep is None or fold is None:
        cv_obj_sub = cv_obj.full
    else:
        cv_obj_sub = cv_obj.cv[rep][fold]
    assert 'eig' in cv_obj_sub
    eig_ev = cv_obj_sub['other']['eig_ev']
    eig_pca_pv = path.join(cv_obj_sub['odir'],'eig.pca.pv')
    cmd = cmd.format(eig_path=eig_path, eig_ev=eig_ev, eig_pca_pv=eig_pca_pv)
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    if p_status == '0':
        eig_k = cv_eig_twstats_k(eig_pca_pv, k_max, alpha)
    return (rep, fold, { 'eig_pca_pv':eig_pca_pv , 'eig_k':eig_k }, cmd, p_stdout, p_status)

"""
"""
#def cv_eig_assoc(cv_obj, eig_path, cmd, k=10, rep=None, fold=None, skip=False):
    #eig_res = eig_log = None
    #p_stdout = 'SKIPPED'; p_status = 0
    #cv_obj_sub = None
    #if rep is None or fold is None:
        #cv_obj_sub = cv_obj.full
    #else:
        #cv_obj_sub = cv_obj.cv[rep][fold]
    #eig_res     = path.join(cv_obj_sub['odir'],'eig.res')
    #eig_log     = path.join(cv_obj_sub['odir'],'eig.log')
    #eig_pca     = cv_obj_sub['other']['eig_pca']
    #eig_geno    = cv_obj_sub['other']['eig_geno']
    #eig_snps    = cv_obj_sub['other']['eig_snps']
    #eig_pheno   = cv_obj_sub['other']['eig_pheno']
    #eig_k       = k #cv_obj_sub['other']['eig_k']
    #cmd = cmd.format( eig_path=eig_path, eig_geno=eig_geno, eig_snps=eig_snps, eig_pheno=eig_pheno, eig_pca=eig_pca, eig_res=eig_res, eig_log=eig_log, eig_k=eig_k)
    #if not skip:
        #p_stdout, p_status = run_cmd(cmd)
    #return (rep, fold, { 'eig_res':eig_res , 'eig_log':eig_log }, cmd, p_stdout, p_status)

"""
"""
#def cv_eig_lambda(cv_obj, eig_path, cmd, rep=None, fold=None, skip=False):
    #eig_lambda = None
    #p_stdout = 'SKIPPED'; p_status = 0
    #cv_obj_sub = None
    #if rep is None or fold is None:
        #cv_obj_sub = cv_obj.full
    #else:
        #cv_obj_sub = cv_obj.cv[rep][fold]
    #eig_res = cv_obj_sub['other']['eig_res']
    #eig_lambda = path.join(cv_obj_sub['odir'],'eig.lambda')
    #cmd = cmd.format(eig_path=eig_path, eig_res=eig_res, eig_lambda=eig_lambda)
    #if not skip:
        #p_stdout, p_status = run_cmd(cmd)
    #return (rep, fold, { 'eig_lambda':eig_lambda }, cmd, p_stdout, p_status)

"""
Run EIGENSTRAT with k <= k_max and compute lambda, stop if lambda close enough to 1
"""
def cv_eig_sel_k(cv_obj, eig_path, eig_assoc_cmd, eig_lambda_cmd, k_max, k_step=1, l_min=1.0, rep=None, fold=None, skip=False):
    info = ''
    cv_obj_sub = None
    p_cmds = []; p_os = []; p_ss = []
    l_p = re.compile('^lambda.*')
    k = 0; l = 100.0
    if rep is None or fold is None:
        cv_obj_sub = cv_obj.full
    else:
        cv_obj_sub = cv_obj.cv[rep][fold]
    assert 'eig' in cv_obj_sub
    eig_sel     = cv_obj_sub['eig']['eig.select']
    eig_lambda  = cv_obj_sub['eig']['eig.lambda']
    eig_res     = cv_obj_sub['eig']['eig.res']
    eig_log     = cv_obj_sub['eig']['eig.log']
    eig_pca     = cv_obj_sub['eig']['eig.pca']
    eig_geno    = cv_obj_sub['eig']['eig.geno']
    eig_snps    = cv_obj_sub['eig']['eig.snps']
    eig_pheno   = cv_obj_sub['eig']['eig.pheno']
    eig_k       = cv_obj_sub['eig']['eig.k']
    
    k_l = {}
    if not skip:
        with open(eig_sel,'w') as res:
            while l >= l_min and l > 1.0 and k < k_max: # stop if l <= 1.0 OR l < l_min OR reached max. number of PCs [l should not be < 1.0 w.r.t. EIGENSTRAT]
                if k == 0: k = 1
                else: k = min(k_max,k+k_step)
                # EIGENSTRAT
                cmd = eig_assoc_cmd.format( eig_path=eig_path, eig_geno=eig_geno, eig_snps=eig_snps, eig_pheno=eig_pheno, eig_pca=eig_pca, eig_res=eig_res, eig_log=eig_log, eig_k=k)
                p_stdout, p_status = run_cmd(cmd)
                p_cmds.append(cmd); p_os.append(p_stdout); p_ss.append(p_status)
                if p_status != 0:  stdout.write('break eig: %s' % p_status); break
                # LAMBDA
                cmd = eig_lambda_cmd.format(eig_path=eig_path, eig_res=eig_res, eig_lambda=eig_lambda)
                p_stdout, p_status = run_cmd(cmd)
                p_cmds.append(cmd); p_os.append(p_stdout); p_ss.append(p_status)
                if p_status != 0:  stdout.write('break lambda: %s' % p_status); break
                # Get LAMBDA
                with open(eig_lambda,'r') as l_file:
                    for line in l_file:
                        if line and l_p.match(line) is not None:
                            #stdout.write("k=%d: %s" %(k,line))
                            l = line.rstrip("\n").split(' ')[1]
                            l = float(l.replace('lambda=',''))
                            break
                res.write("%d\t%.5f\n" % (k,l)); res.flush()
                k_l[k] = l
        # Reached max PC num. but lambda still "bad" -> get k with min. lambda
        if k == k_max and l >= l_min and l > 1.0:
            k = 1; l = k_l[k]
            for k_ in sorted(k_l.keys()):
                if k_l[k_] < k_l[k]: k = k_
            info += "Minimal k with minimal lambda is %d with %.5f\n" %(k, k_l[k])
            # EIGENSTRAT
            cmd = eig_assoc_cmd.format( eig_path=eig_path, eig_geno=eig_geno, eig_snps=eig_snps, eig_pheno=eig_pheno, eig_pca=eig_pca, eig_res=eig_res, eig_log=eig_log, eig_k=k)
            p_stdout, p_status = run_cmd(cmd)
            p_cmds.append(cmd); p_os.append(p_stdout); p_ss.append(p_status)
            if p_status != 0:  stdout.write('break eig: %s' % p_status)
            # LAMBDA
            cmd = eig_lambda_cmd.format(eig_path=eig_path, eig_res=eig_res, eig_lambda=eig_lambda)
            p_stdout, p_status = run_cmd(cmd)
            p_cmds.append(cmd); p_os.append(p_stdout); p_ss.append(p_status)
            if p_status != 0:  stdout.write('break lambda: %s' % p_status)
        info += "\n".join(['%s : %s\n%s\n' % (p_cmds[i],p_ss[i],p_os[i]) for i in range(0,len(p_cmds))])
        with open(eig_k,'w') as output: output.write("%d\n" % k)
    else: info = '\nSKIPPED\n'
    return (rep, fold, info)

#def eig_assoc_lambda(cmd_assoc, cmd_lambda, eig_lambda, k):
    #reurn "Removed"
    #l_p = re.compile('^lambda.*')
    #l = None
    ## ASSOC
    #p_stdout, p_status = run_cmd(cmd_assoc)
    #if p_status != 0:  stdout.write('break eig: %s:\n %s' % (cmd_assoc, p_stdout)); stdout.flush()
    ## LAMBDA
    #p_stdout2, p_status2 = run_cmd(cmd_lambda)
    #if p_status2 != 0:  stdout.write('break lambda: %s:\n %s' % (cmd_lambda, p_stdout2)); stdout.flush()
    ## Get LAMBDA
    #stdout.write(cmd_lambda); stdout.flush()
    #with open(eig_lambda,'r') as l_file:
        #for line in l_file:
            #if line and l_p.match(line) is not None:
                ##stdout.write("k=%d: %s" %(k,line))
                #l = line.rstrip("\n").split(' ')[1]
                #l = float(l.replace('lambda=',''))
                #break
    #return (cmd_assoc, p_stdout, p_status, cmd_lambda, p_stdout2, p_status2, k, l)
"""
"""
#def cv_eig_sel_k_parallel(cv_obj, eig_path, eig_assoc_cmd, eig_lambda_cmd, k_max, k_step=1, l_min=1.0, rep=None, fold=None, cores=1, skip=False):
    #reurn "Removed"
    #if skip: return "Sel. k PCs w.r.t. lambda: SKIPPED\n"
    #cv_obj_sub = None
    #if rep is None or fold is None:
        #cv_obj_sub = cv_obj.full
    #else:
        #cv_obj_sub = cv_obj.cv[rep][fold]
    #assert 'eig' in cv_obj_sub
    #eig_sel     = cv_obj_sub['eig']['eig.select']
    #eig_lambda  = cv_obj_sub['eig']['eig.lambda']
    #eig_res     = cv_obj_sub['eig']['eig.res']
    #eig_log     = cv_obj_sub['eig']['eig.log']
    #eig_pca     = cv_obj_sub['eig']['eig.pca']
    #eig_geno    = cv_obj_sub['eig']['eig.geno']
    #eig_snps    = cv_obj_sub['eig']['eig.snps']
    #eig_pheno   = cv_obj_sub['eig']['eig.pheno']
    #eig_k       = cv_obj_sub['eig']['eig.k']
    
    #tmp_dir = path.join(cv_obj_sub['odir'],'tmp')
    #if not path.exists(tmp_dir): os.makedirs(tmp_dir)
    #eig_res_k        = [ path.join(tmp_dir,"eig.res.%d" % k)   for k in range(1,k_max+1)]
    #eig_log_k        = [ path.join(tmp_dir,"eig.log.%d" % k)   for k in range(1,k_max+1)]
    #eig_lambda_k     = [ path.join(tmp_dir,"eig.lamda.%d" % k) for k in range(1,k_max+1)]
    #eig_assoc_cmd_k  = [None] * len(range(1,k_max+1))
    #eig_lambda_cmd_k = [None] * len(range(1,k_max+1))
    #for k in range(1,k_max+1):
        #eig_assoc_cmd_k[k-1]  = eig_assoc_cmd.format( eig_path=eig_path, eig_geno=eig_geno, eig_snps=eig_snps, eig_pheno=eig_pheno, eig_pca=eig_pca, eig_res=eig_res_k[k-1], eig_log=eig_log_k[k-1], eig_k=k)
        #eig_lambda_cmd_k[k-1] = eig_lambda_cmd.format(eig_path=eig_path, eig_res=eig_res_k[k-1], eig_lambda=eig_lambda_k[k-1])
        
    #pool        = Pool(cores)
    #pool_iter   = itertools.product(eig_assoc_cmd_k,eig_lambda_cmd_k,eig_lambda_k,range(1,k_max+1))
    #results     = pool.starmap( eig_assoc_lambda , pool_iter )
    #pool.close(); pool.join()
    
    #k_l = {}
    #for res in results:
        #k_l[res[6]] = res[7] # save k : l
        ## add info for log
        #info += '%s : %s\n%s\n' % (res[0],res[2],res[1])
        #info += '%s : %s\n%s\n' % (res[3],res[5],res[4])
    ## remove files
    #if path.isdir(tmp_dir): shutil.rmtree(tmp_dir)
    ## select k
    #sel_k = 1
    #with open(eig_sel,'w') as res:
        #for k in sorted(k_l.keys()):
            #res.write("%d\t%.5f\n" % (k,k_l[k])); res.flush()
            #if k_l[k] == 1.0 or k_l[k] < l_min: # condition is fullfilled
                #sel_k = k
                #break
            #elif k_l[k] < k_l[sel_k]: # found better lambda
                #sel_k = k
    #info += "Selected k=%d with lambda = %.5f\n" % (sel_k, k_l[sel_k])
    #with open(eig_k,'w') as output: output.write("%d\n" % k)
    ## run ASSOC and LAMBDA witk selected k
    ## EIGENSTRAT
    #cmd = eig_assoc_cmd.format( eig_path=eig_path, eig_geno=eig_geno, eig_snps=eig_snps, eig_pheno=eig_pheno, eig_pca=eig_pca, eig_res=eig_res, eig_log=eig_log, eig_k=sel_k)
    #p_stdout, p_status = run_cmd(cmd)
    #info += '%s : %s\n%s\n' % (cmd,p_status,p_stdout)
    #if p_status != 0:  stdout.write('break eig: %s' % p_status)
    ## LAMBDA
    #cmd = eig_lambda_cmd.format(eig_path=eig_path, eig_res=eig_res, eig_lambda=eig_lambda)
    #p_stdout, p_status = run_cmd(cmd)
    #info += '%s : %s\n%s\n' % (cmd,p_status,p_stdout)
    #if p_status != 0:  stdout.write('break lambda: %s' % p_status)
    #return (rep, fold, info)

def cv_eig_post_proc(cv_obj, src_path, cmd, rep, fold, skip=False):
    assert 'eig' in cv_obj.cv[rep][fold]
    p_stdout = 'SKIPPED'; p_status = 0
    
    eig_res     = cv_obj.cv[rep][fold]['eig']['eig.res']
    eig_snps    = cv_obj.cv[rep][fold]['eig']['eig.snps']
    eig_ev      = cv_obj.cv[rep][fold]['eig']['eig.ev']
    eig_k       = cv_obj.cv[rep][fold]['eig']['eig.k']
    eig_sel     = cv_obj.cv[rep][fold]['features_sel']
    
    k = None
    with open(eig_k,'r') as input:
        line = input.readline()
        line = line.rstrip("\n")
        k = int(line)
    
    cmd = cmd.format( src_path=src_path, eig_res=eig_res, eig_snps=eig_snps, eig_ev=eig_ev, eig_k=k, features_sel=eig_sel )
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (rep, fold, cmd, p_stdout, p_status)

def cv_eig_full_post_proc(cv_obj, src_path, cmd, cores=2, skip=False):
    assert 'eig' in cv_obj.full
    p_stdout = 'SKIPPED'; p_status = 0
    
    eig_res     = cv_obj.full['eig']['eig.res']
    eig_snps    = cv_obj.full['eig']['eig.snps']
    eig_ev      = cv_obj.full['eig']['eig.ev']
    eig_k       = cv_obj.full['eig']['eig.k']

    k = None
    with open(eig_k,'r') as input:
        line = input.readline()
        line = line.rstrip("\n")
        k = int(line)
    
    cmd = cmd.format( \
        src_path=src_path, \
        eig_res=eig_res, eig_snps=eig_snps, eig_ev=eig_ev, eig_k=k, \
        geno_file=cv_obj.X_bin, pheno_file=cv_obj.Y_file, pheno_name=cv_obj.Y_name, \
        samples=cv_obj.full['samples'], \
        cores=cores \
    )
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (cmd, p_stdout, p_status)

"""
Check sample lists
"""
def cv_check_samples(in_samples_f, ex_samples_f, cv_obj):
    if cv_obj.done['sample_check']:
        return "SKIP: Samples already checked"
    info = ""
    in_samples = set(read_id_list(in_samples_f))
    ex_samples = set(read_id_list(ex_samples_f))
    
    # CV full:
    cv_samples = set(read_id_list(cv_obj.full['samples']))
    # no intersectoion with excluded
    assert len(cv_samples.intersection(ex_samples))==0, "Intersection with excluded: Assertion error in %s" % cv_samples
    # all in included
    assert len(cv_samples.difference(in_samples))  ==0, "Set diff. with included: Assertion error in %s" % cv_samples
    info += "CV %s (full): Sample list %s: %d entries\n" % (cv_obj.Y_name,cv_obj.full['samples'],len(cv_samples))
    
    # CV folds: only samples from samples in CV obj., print number
    for rep in cv_obj.cv.keys():
        for fold in cv_obj.cv[rep].keys():
            train_samples = set(read_id_list(cv_obj.cv[rep][fold]['samples_train']))
            test_samples  = set(read_id_list(cv_obj.cv[rep][fold]['samples_test']))
            # no intersection of train and test
            assert len(train_samples.intersection(test_samples))==0, "Assertion error in %s and %s" % (train_samples, test_samples)
            # together they should be the same as in full set
            assert cv_samples == train_samples.union(test_samples), "Assertion error in %s and %s" % (train_samples, test_samples)
            # number of samples
            info += "CV %s (rep %d, fold %d): Sample list %s: %d entries\n" % (cv_obj.Y_name,rep,fold,cv_obj.cv[rep][fold]['samples_train'],len(train_samples))
            info += "CV %s (rep %d, fold %d): Sample list %s: %d entries\n" % (cv_obj.Y_name,rep,fold,cv_obj.cv[rep][fold]['samples_test'],len(test_samples))
    cv_obj.done['sample_check'] = True
    return info

def cv_model(cv_obj, src_path, model, model_cmd, model_params, rep=None, fold=None, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    if rep is None or fold is None:
        x_file = cv_obj.X_bin
        y_file = cv_obj.Y_file
        y_pheno= cv_obj.Y_name
        o_dir  = cv_obj.full['odir']
        o_bname= model
        samples_train = cv_obj.full['samples']
        features      = cv_obj.full['features_sel']
        cmd = model_cmd.format( \
            src_path=src_path, \
            x_file=x_file, y_file=y_file, y_pheno=y_pheno, \
            o_dir=o_dir, o_bname=o_bname, \
            samples_train=samples_train, features=features, \
            model_params=model_params \
        )
    else:
        x_file = cv_obj.X_bin
        y_file = cv_obj.Y_file
        y_pheno= cv_obj.Y_name
        o_dir  = cv_obj.cv[rep][fold]['odir']
        o_bname= model
        samples_train = cv_obj.cv[rep][fold]['samples_train']
        samples_test  = cv_obj.cv[rep][fold]['samples_test']
        features      = cv_obj.cv[rep][fold]['features_sel']
        cmd = model_cmd.format( \
            src_path=src_path, \
            x_file=x_file, y_file=y_file, y_pheno=y_pheno, \
            o_dir=o_dir, o_bname=o_bname, \
            samples_train=samples_train, samples_test=samples_test, features=features, \
            model_params=model_params \
        )
    if not skip:
        ## delay:
        #if rep is not None or fold is not None:
            #time.sleep((fold-1)*3 + (rep-1)*3 + 5)
        p_stdout, p_status = run_cmd(cmd)
    return (rep, fold, cmd, p_stdout, p_status)

"""
"""
def cv_fold_sum(cv_obj, src_path, model, rep, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    
    pred_files =  [ fold[model]["pred.csv"] for fold in cv_obj.cv[rep].values() ]
    cmd = "Rscript {src_path}/utils/cv_fold_sum.R -pred_files {pred_files} -o_file {o_file} --plot --src_path {src_path}/utils --verbose"
    
    cmd = cmd.format(src_path=src_path, pred_files=' '.join(pred_files), o_file=path.join(cv_obj.odir,"%s_rep%d_perf.csv" % (model,rep)))
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (rep, cmd, p_stdout, p_status)

"""
"""
def cv_rep_sum(cv_obj, src_path, model, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    
    i_files =  [ path.join(cv_obj.odir,"%s_rep%d_perf.csv" % (model,rep)) for rep in cv_obj.cv.keys() ]
    cmd = "Rscript {src_path}/utils/cv_rep_sum.R -perf_files {i_files} -o_file {o_file} --src_path {src_path}/utils --reps {reps} --min_rep_pct 50.0 --verbose"
    
    cmd = cmd.format(src_path=src_path, i_files=' '.join(i_files), o_file=path.join(cv_obj.odir,"%s_total_perf.csv" % model), reps=cv_obj.reps)
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (cmd, p_stdout, p_status)

"""
"""
def cv_pheno_sum(cv_objs, src_path, odir, model, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    
    i_perf  = [ path.join(cv_obj.odir,"%s_total_perf.csv" % model) for cv_obj in cv_objs if cv_obj.Y_check ]
    i_names = [ cv_obj.Y_name for cv_obj in cv_objs if cv_obj.Y_check ]
    i_mods  = [ cv_obj.full[model]['combi.txt'] for cv_obj in cv_objs if cv_obj.Y_check ]
    o_file  = path.join(odir,'%s_CV_summary.csv' % model)
    
    cmd = "Rscript {src_path}/utils/cv_y_sum.R -y_perf_files {i_perf} -y_names {i_names} -y_models {i_mods} -o_file {o_file} --verbose"
    
    cmd = cmd.format(src_path=src_path, i_perf=' '.join(i_perf), i_names=' '.join(i_names), i_mods=' '.join(i_mods), o_file=o_file)
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (cmd, p_stdout, p_status)

"""
"""
def cv_combi_sel_features(cv_obj, src_path, skip=False):
    p_stdout = 'SKIPPED'; p_status = 0
    
    i_files = []
    for rep in cv_obj.cv.keys():
        for fold in cv_obj.cv[rep].keys():
            i_files.append(cv_obj.cv[rep][fold]['features_sel'])
    
    cmd = "Rscript {src_path}/utils/combine_sel_features.R -sel_features_files {i_files} -o_file {o_file} --src_path {src_path}/utils --verbose"
    cmd = cmd.format(src_path=src_path, i_files=' '.join(i_files), o_file=cv_obj.full['features_sel'])
    
    if not skip:
        p_stdout, p_status = run_cmd(cmd)
    return (cmd, p_stdout, p_status)

def cv_clean(cv_obj):
    info = "Removing:\n"
    # full
    for f_key, f in cv_obj.full['eig'].items():
        if isinstance(f, str) and path.isfile(f) and any([p.match(f) is not None for p in rm_files]):
            info += "\t%s\n" % f
            if path.isfile(f): os.remove(f)
    # folds
    for rep in cv_obj.cv.keys():
        for fold in cv_obj.cv[rep].keys():
            for f_key, f in cv_obj.cv[rep][fold]['eig'].items():
                if isinstance(f, str) and path.isfile(f) and any([p.match(f) is not None for p in rm_files]):
                    info += "\t%s\n" % f
                    if path.isfile(f): os.remove(f)
    # NOTE hard-coded
    cv_obj.done['full_eig_conv'] = False
    cv_obj.done['cv_eig_conv']   = False
    return info
