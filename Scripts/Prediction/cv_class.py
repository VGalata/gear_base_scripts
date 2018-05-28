#!/bin/python
# -*- coding: utf-8 -*-

import sys
from os import path
import itertools
from multiprocessing import Pool
import shlex, subprocess
import copy

from utils.my_utils import make_str_compliant

class CVClass:
    """Prediction class to store relevant information as feature and phenotype data, output files etc."""
    
    def __init__(self, X_file, X_type, X_bin, X_gds, Y_file, Y_type, Y_name, odir, X_vcf_add=None, in_samples=None, ex_samples=None, features=None, cv_reps=1, cv_folds=1):
        # Features
        self.X_file     = X_file    # given feature file
        self.X_type     = X_type    # type of feature file
        self.X_bin      = X_bin     # feature file in binary format (e.g. if VCF converted to GDS and than saved as binary matrix)
        self.X_gds      = X_gds     # VCF converted to GDS
        self.X_vcf_add  = X_vcf_add # additional file if X is a VCF (gene presence/absence matrix needed later for feature pre-processing)
        
        # Phenotype
        self.Y_file     = Y_file    # phenotype file
        self.Y_type     = Y_type    # phenotype file format/type
        self.Y_name     = Y_name    # phenotype name, i.e. column name in the phenotype file
        self.Y_name_str = make_str_compliant(Y_name) # name made compliant to be used as dir./file name
        self.Y_stat     = dict.fromkeys(['file','total','miss','class_num','class_pct'], None) # phenotype stat.s, e.g. class size, total, etc.
        self.Y_check    = False # True if enough samples and more or less balanced classes (bin. Y)
        
        # Output dir.
        self.odir       = path.join(odir,self.Y_name_str) # Output dir.
        
        # Samples/features
        self.in_samples     = in_samples # given sample lists
        self.ex_samples     = ex_samples # given sample lists
        self.samples        = path.join(self.odir,'samples.txt') # samples in X and Y with non-miss. value in Y
        self.Y_stat['file'] = path.join(self.odir,'pheno_stats.csv')
        self.features       = features # all features
        
        # Full instance
        self.full = {}
        self.full['bname']          = 'full' # basename, kind of a run ID <-> CV
        self.full['odir']           = path.join(self.odir,self.full['bname'])
        self.full['samples']        = self.samples
        self.full['features']       = self.features
        self.full['features_pr']    = path.join(self.full['odir'],'features_pr.txt')
        self.full['features_sel']   = path.join(self.full['odir'],'features_sel.txt')
        
        # k-fold CV instances for each rep.
        self.reps = cv_reps
        self.folds= cv_folds
        self.cv = {}
        for rep in range(1,cv_reps+1):
            rep_dict = {}
            for fold in range(1,cv_folds+1):
                fold_dict = {}
                fold_dict['bname']          = 'rep%d_fold%d' % (rep, fold)
                fold_dict['odir']           = path.join(self.odir,fold_dict['bname'])
                fold_dict['samples_train']  = path.join(fold_dict['odir'],'samples_train.txt') # training samples
                fold_dict['samples_test']   = path.join(fold_dict['odir'],'samples_test.txt') # testing samples
                fold_dict['features']       = self.features
                fold_dict['features_pr']    = path.join(fold_dict['odir'],'features_pr.txt')
                fold_dict['features_sel']   = path.join(fold_dict['odir'],'features_sel.txt')
                rep_dict[fold] = fold_dict
            self.cv[rep] = rep_dict
        
        # Dictionary to know which steps were already performed
        # TODO: steps regarding the model: should depend on model name
        self.done = {\
            'samples':False, 'sample_check':False,\
            'y_stats':False, 'y_check':False,\
            'cv_folds':False,\
            'cv_preproc':False,\
            'cv_eig_conv':False, 'cv_eig_pca':False, 'cv_eig_assoc':False, 'cv_eig_proc':False,\
            'full_preproc':False,\
            'full_eig_conv':False, 'full_eig_pca':False, 'full_eig_assoc':False, 'full_eig_proc':False, 'full_sel_f':False\
        }
    
    # NOTE dict.s for EIGENSTRAT/models
    eig_dict   = dict.fromkeys( ['eig.geno','eig.snps','eig.snps.rm','eig.pheno','eig.pca','eig.pca.evec','eig.pca.log','eig.pca.par','eig.pca.plot','eig.ev','eig.res','eig.res.adj','eig.res.par','eig.log','eig.select','eig.k','eig.lambda'] )
    mod_dict   = dict.fromkeys( ['rds','pred.csv','combi.txt','pdf'] )
    mod_dict_done = {'cv_mod':False, 'cv_fold_sum':False, 'cv_rep_sum':False, 'full_mod':False}
    
    def init_dict(self, d, odir):
        for k in d.keys():
            d[k] = path.join(odir,k)
        return d
    
    def init_mod_dict(self, d, odir, mod_name):
        for k in d.keys():
            d[k] = path.join(odir,"%s.%s" % (mod_name,k))
        return d
    
    def update_dict(self, d_old, d_init):
        d_init.update(d_old) # overwrite by old values
        d_old.update(d_init) # will add attr. contained in init. dict. but not in old
        return d_old
    
    def add_eig(self):
        if 'eig' not in self.full:
            self.full['eig'] = self.init_dict(copy.deepcopy(self.eig_dict),self.full['odir'])
        for rep in self.cv.keys():
            for fold in self.cv[rep].keys():
                if 'eig' not in self.cv[rep][fold]:
                    self.cv[rep][fold]['eig'] = self.init_dict(copy.deepcopy(self.eig_dict),self.cv[rep][fold]['odir'])
    
    def add_mod(self, mod_name):
        if mod_name not in self.full:
            self.full[mod_name] = self.init_mod_dict(copy.deepcopy(self.mod_dict),self.full['odir'], mod_name)
        for rep in self.cv.keys():
            for fold in self.cv[rep].keys():
                if mod_name not in self.cv[rep][fold]:
                    self.cv[rep][fold][mod_name] = self.init_mod_dict(copy.deepcopy(self.mod_dict),self.cv[rep][fold]['odir'], mod_name)
        if mod_name not in self.done:
            self.done[mod_name] = copy.deepcopy(self.mod_dict_done)
    
    # NOTE ???
    def set_cv_attr_to_value(self, attr, value):
        for rep in self.cv.keys():
            for fold in self.cv[rep].keys():
                self.cv[rep][fold][attr] = value
    
    # NOTE Print methods
    def print_instance(self, inst):
        s = "Instance %s: %s\n" % (inst['bname'],inst['odir'])
        attr_print = set(inst.keys())
        attr_print.difference_update(['odir','bname','other'])
        for attr in attr_print:
            s += "\t%s: %s\n" % (attr,inst[attr])
        return s
    
    def print_class_num(self):
        if self.Y_stat['class_num'] is None:
            return ''
        else:
            return ' ; '.join([' - '.join([k,"%d" % v]) for k,v in self.Y_stat['class_num'].items()])
    
    def print_class_pct(self):
        if self.Y_stat['class_pct'] is None:
            return ''
        else:
            return ' ; '.join([' - '.join([k,"%.2f" % v]) for k,v in self.Y_stat['class_pct'].items()])
    
    def print_class_stat(self):
        s = "Y stats: Total="
        if self.Y_stat['total'] is None:
            s += "<None>"
        else:
            s += "%d" % self.Y_stat['total']
        s += ", Miss="
        if self.Y_stat['miss'] is None:
            s += "<None>"
        else:
            s += "%d" % self.Y_stat['miss']
        s += ", Class count: %s, Class pct: %s\n"  % (self.print_class_num(), self.print_class_pct())
        return s
    
    def print_cv_class(self):
        s  = "#%s#\n" % ("-"*50)
        s += "CV object %s:\n" % self.Y_name
        s += "X file: %s\n\tType: %s\n\tBin. file: %s\n\tGDS file: %s\n"  % (self.X_file, self.X_type, self.X_bin, self.X_gds)
        s += "Y file: %s\n\tType: %s\n\tName: %s\n\tOutput name: %s\n"  % (self.Y_file, self.Y_type, self.Y_name, self.Y_name_str)
        s += self.print_class_stat()
        s += "Samples to include: %s\n" % (self.in_samples)
        s += "Samples to exclude: %s\n" % (self.ex_samples)
        s += "Samples: %s\n" % self.samples
        s += "%d-fold CV with %d repetitions\n" % (self.folds,self.reps)
        s += self.print_instance(self.full)
        #for rep in self.cv.keys():
            #for fold in self.cv[rep].keys():
                #s += self.print_instance(self.cv[rep][fold])
        s  += "#%s#\n" % ("-"*50)
        return s
