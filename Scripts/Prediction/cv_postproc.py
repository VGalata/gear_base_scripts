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
import glob
import pickle
#import shutil
#import logging
import re
#import time, datetime
import copy
import argparse
import itertools
from multiprocessing import Pool
import shlex, subprocess

from cv_class import CVClass
from utils.my_utils import read_id_list
from utils.bio_utils import annot_from_gff

gene_annot_fields = ['gene','product','gene_length','strand']
vcf_annot_fields  = ['gene','product','gene_length','strand','chrom','pos','ref','alt','geneID','DNA_change','AA_change']

def load_cv_obj(cv_file, verbose=True):
    with open(cv_file, 'rb') as input:
        return( pickle.load(input) )

def read_lambda(eig_lambda):
    l_p = re.compile('^lambda.*')
    with open(eig_lambda,'r') as l_file:
        for line in l_file:
            if line and l_p.match(line) is not None:
                l = line.rstrip('\n').split(' ')[1]
                return ( float(l.replace('lambda=','')) )

def read_k(eig_k):
    with open(eig_k,'r') as input:
        line = input.readline()
        line = line.rstrip('\n')
        return ( int(line) )

def pheno_general_stats(obj, model):
    info = ""
    obj_ = load_cv_obj(obj)
    info += "%s\nPehnotype: %s\n" % (sep_str,obj_.Y_name)
    # number of samples (after sample intersection)
    info += "\tNumber of samples without missing pheno: %d\n" % len(set(read_id_list(obj_.samples)))
    # pheno. stats (after sample intersection)
    info += "\t0:%d,1:%d | 0:%.2f,1:%.2f\n" % (obj_.Y_stat['class_num']['0'],obj_.Y_stat['class_num']['1'],obj_.Y_stat['class_pct']['0'],obj_.Y_stat['class_pct']['1'])
    # pheno check
    if obj_.Y_check:
        # k and lambda
        k = read_k(obj_.full['eig']['eig.k'])
        l = read_lambda(obj_.full['eig']['eig.lambda'])
        info += "EIGENSTRAT\n\tFull: k=%d, lambda=%.5f\n" % (k,l)
        for rep in obj_.cv.keys():
            info += "\tRep %d" % (rep)
            for fold in obj_.cv[rep].keys():
                k = read_k(obj_.cv[rep][fold]['eig']['eig.k'])
                l = read_lambda(obj_.cv[rep][fold]['eig']['eig.lambda'])
                info += " | fold %d: k=%d, lambda=%.5f" % (fold,k,l)
            info += "\n"
        # CV selected features
        n = len(set(read_id_list(obj_.full['features_sel'])))
        info += "Selected features:\n\tFull: %d\n" % n
        for rep in obj_.cv.keys():
            info += "\tRep %d" % (rep)
            for fold in obj_.cv[rep].keys():
                n = len(set(read_id_list(obj_.cv[rep][fold]['features_sel'])))
                info += " | fold %d: %d" % (fold,n)
            info += "\n"
        # CV model
        s = ';'.join(list(set(read_id_list(obj_.full[model]["combi.txt"]))))
        info += "Models\n\tFull: %s\n" % s
        for rep in obj_.cv.keys():
            for fold in obj_.cv[rep].keys():
                s = ';'.join(list(set(read_id_list(obj_.cv[rep][fold][model]["combi.txt"]))))
                info += "\tRep %d, fold %d: %s\n" % (rep,fold,s)
        # CV eval. reps
        info += 'Mssing fold pred. results:'
        for rep in obj_.cv.keys():
            for fold in obj_.cv[rep].keys():
                f = obj_.cv[rep][fold][model]['pred.csv']
                if not os.path.isfile(f) or os.stat(f).st_size == 0:
                    info += "Rep %d, fold %d: no model prediction results\n" % (rep,fold)
        info += '\n'
    else:
        info += 'Check not passed\n'
    info += "%s\n" % sep_str
    return info

#def read_gff_annot_tab(gff_annot_tab):
    #annot = {}
    #annot_fields = None
    #header = True
    #with open(gff_annot_tab,'r') as gff_annot:
        #for line in gff_annot:
            #line = line.rstrip('\n')
            #if not line: continue
            #if header:
                #header = False
                #annot_fields = line.split('\t')[1:] # fist column contains ID
                #continue
            #line = line.split('\t') # split info fields
            #assert line[0] not in annot
            #assert annot_fields is not None
            #annot[line[0]] = dict.fromkeys(annot_fields,None) # create empty dir.
            #for i in range(0,len(line)):
                #annot[line[0]][annot_fields[i-1]] = line[i]
    #return annot

def read_annot_tab(annot_tab, id_prefix=None):
    annot = {}
    annot_fields = None
    header = True
    with open(annot_tab,'r') as tab_annot:
        for line in tab_annot:
            line = line.rstrip('\n')
            if not line: continue
            if header:
                header = False
                annot_fields = line.split('\t')[1:] # fist column contains ID
                continue
            line = line.split('\t') # split info fields
            if id_prefix is not None: # add id prefix if given
                line[0] = "%s%s" % (id_prefix,line[0])
            assert line[0] not in annot
            assert annot_fields is not None
            annot[line[0]] = dict.fromkeys(annot_fields,None) # create empty dir.
            for i in range(1,len(line)):
                annot[line[0]][annot_fields[i-1]] = line[i]
    return annot
    

def add_gwas_annot(obj, annot, odir, obname, model, f_type):
    fields = None
    if f_type == 'gff':
        fields = gene_annot_fields
    elif f_type == 'vcf':
        fields = vcf_annot_fields
    obj_ = load_cv_obj(obj)
    if not obj_.Y_check: return
    
    gwas_res = obj_.full['eig']['eig.res.adj']
    if not path.isfile(gwas_res) or os.stat(gwas_res).st_size == 0: return
    
    ofile = path.join(odir,"%s_%s_gwas.csv" % (obname, obj_.Y_name_str))
    ID_col = 2
    with open(gwas_res,'r') as gwas_i, open(ofile,'w') as gwas_o:
        header = True
        for line in gwas_i:
            line   = line.rstrip('\n')
            line_s = line.split('\t')
            if not line: continue
            if header:
                header = False
                gwas_o.write("%s\t%s\n" % (line,'\t'.join(fields)))
                if line_s[2] != 'Feature':
                    ID_col = [i for i in range(0,len(line_s)) if line_s[i]=='Feature'][0]
                continue
            ID = line_s[ID_col]
            assert ID in annot, "Feature %s not in created annotation dictionary" % ID
            gwas_o.write("%s\t%s\n" % (line,'\t'.join([ annot[ID][f] for f in fields ])))

def add_model_annot(obj, annot, odir, obname, model, f_type):
    fields = None
    if f_type == 'gff':
        fields = gene_annot_fields
    elif f_type == 'vcf':
        fields = vcf_annot_fields
    
    obj_ = load_cv_obj(obj)
    if not obj_.Y_check: return
    
    model_res = path.join(obj_.odir,"%s_total_perf.csv" % model)
    if not path.isfile(model_res) or os.stat(model_res).st_size == 0: return # no CV perf.
    if not path.isfile(obj_.full[model]['combi.txt']): return # no/empty model
    
    ofile = path.join(odir,"%s_%s_%s.csv" % (obname, obj_.Y_name_str, model))
    perf_list = ['ERR','ACC','B_ACC','SENS','SPEC','PREC','NPV','FPR','FNR','Fmeasure','gm_RS','gm_RP','AUC_ROC','AUC_PR']
    perf_dict = dict.fromkeys( perf_list , None )
    # read perf
    with open(model_res) as ifile:
        header = True
        for line in ifile:
            if not line: continue
            line = line.rstrip('\n')
            line = line.split('\t')
            if header:
                header = False
                for i in range(0,len(line)):
                    if line[i] in perf_dict:
                        if perf_dict[line[i]] is None: perf_dict[line[i]] = {}
                        perf_dict[line[i]]['id'] = i
                continue
            if line[0] == 'mean':
                for k in perf_dict.keys():
                    perf_dict[k]['value'] = line[perf_dict[k]['id']]
    # read features
    features = set( read_id_list(obj_.full[model]['combi.txt']) )
    # write file
    with open(ofile, 'w') as of:
        of.write("Feature\t%s\t%s\n" % ('\t'.join(fields),'\t'.join(perf_list)))
        for feature in features:
            of.write("%s\t%s\t%s\n" % (feature, '\t'.join([ annot[feature][f] for f in fields ]), '\t'.join([ perf_dict[k]['value'] for k in perf_list ])))

# NOTE MAIN
if __name__ == "__main__":
    #--------------------------------------------------#
    # NOTE Args
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--idir',      '-i', help='', required=True)
    #parser.add_argument('--odir',       '-o', help='', required=True)
    parser.add_argument('--obname',     '-n', help='', required=True)
    parser.add_argument('--model',      '-m', help='', required=True)
    parser.add_argument('--f_type',     '-f', help='', choices=['gff','vcf'], required=True)
    parser.add_argument('--f_source',   '-s', help='GFF annot table or VCF file', required=True)
    parser.add_argument('--cores',      '-c', help='', default=1, type=int)
    parser.add_argument("--verbose",    '-v', help="Verbosity", action="store_true")
    args = parser.parse_args()
    #--------------------------------------------------#
    # NOTE Output dir.
    odir = path.join(args.idir,'Final')
    if not os.path.exists(odir): os.makedirs(odir)
    #--------------------------------------------------#
    # NOTE All pheno. CV objects
    objs = glob.glob("%s/obj/*.pkl" % args.idir)
    #--------------------------------------------------#
    # NOTE Collect information of each object
    sep_str = "#%s#" % ('-'*50)
    with open(path.join(odir,"%s_stats.txt" % args.obname),'w') as ofile:
        for obj in objs:
            info = pheno_general_stats(obj, args.model)
            ofile.write(info)
    #--------------------------------------------------#
    # NOTE Feature annotation
    annot = None
    # All genes
    obj_ = load_cv_obj(objs[0]) # load only first object (all should have same feature list)
    features = set(read_id_list(obj_.features))
    sys.stdout.write("There are %d features\n" % len(features))
    # Annotation from table
    if args.f_type == 'gff':
        annot = read_annot_tab(annot_tab=args.f_source)
    elif args.f_type == 'vcf':
        annot = read_annot_tab(annot_tab=args.f_source, id_prefix='var_')
    assert len(annot) == len(features)
    sys.stdout.write('Annotations were collected\n')
    #--------------------------------------------------#
    # NOTE Add annotation to results
    assert annot is not None
    pool = Pool(args.cores)
    pool_iter = itertools.product(objs,[annot],[odir],[args.obname],[args.model],[args.f_type])
    annots = pool.starmap( add_gwas_annot , pool_iter )
    pool.close(); pool.join()
    # Model results
    pool = Pool(args.cores)
    pool_iter = itertools.product(objs,[annot],[odir],[args.obname],[args.model],[args.f_type])
    annots = pool.starmap( add_model_annot , pool_iter )
    pool.close(); pool.join()
    #--------------------------------------------------#
    