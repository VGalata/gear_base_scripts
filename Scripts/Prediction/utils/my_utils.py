#!/bin/python
# -*- coding: utf-8 -*-

# NOTE Imports
import sys
import os
from os import path
import argparse
import time, datetime
import re
from collections import OrderedDict
import unicodedata
import shlex, subprocess

"""
Get timestamp for printing (YYYY.MM.DD HH:MM:SS) or to add as part of file name (YYYY-MM-DD_HH_MM)
"""
def timestamp():
    return str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d %H:%M:%S'))

def timestamp_():
    return str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M'))

"""
"""
def run_cmd(cmd):
    p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p_stdout = p.stdout.read().decode()
    p_comm   = p.communicate()[0]
    p_status = p.returncode
    return p_stdout, p_status

"""
Make given string compliant for using it as dir./file name
"""
def make_str_compliant(s):
    # UTF8 -> ASCII; replace bad chars
    # http://python-gtk-3-tutorial.readthedocs.org/en/latest/unicode.html
    #return re.sub('/| ','_',unicodedata.normalize('NFKD', s.decode('utf8')).encode('ascii','ignore'))
    return re.sub('/| ','_',unicodedata.normalize('NFKD', s).encode('ascii','ignore').decode())

"""
Print long arrays/lists
"""
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def print_list(l, n, s1='\n', s2=', '):
    if type(l) is not list:
        l = list(l)
    return s1.join([ s2.join(sublist) for sublist in chunks(l,n) ])

"""
Custom parser: Can parse config files (ignores comment lines)
"""
class CustomParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        #if arg_line.find("#",0,1):#ignore comment lines
            #print "hier"
            #yield []
        for arg in arg_line.split():
            #print arg
            if arg_line.find("#",0,1) != -1:#ignore comment lines
                continue
            if not arg.strip():
                continue
            if arg[0] == "\"":
                #combined = arg_line[arg_line.find("\"") +1 :arg_line.rfind("\"")]
                ##print combined
                ##combined = combined.replace("\"", "")
                #yield combined
                #break#don't parse remaining stuff as single args in the current line
                combined = arg_line[arg_line.find("\""):arg_line.rfind("\"")+1]
                for comb_arg in combined.split('\" \"'):
                    yield comb_arg.replace('\"','')
                break
            print(arg)
            yield arg
    
    def print_args(self,args):
        s="Arguments:\n"
        for arg, arg_v in vars(args).items():
            s += "%s: %s\n" % (arg, str(arg_v))
        s += '\n'
        return s

"""
Read list of IDs (one per row)
"""
def read_id_list(list_file):
    l = []
    if list_file is not None and path.isfile(list_file):
        with open(list_file, "r") as list_f:
            for line in list_f:
                line = line.rstrip('\n')
                if line: l.append(line)
    return l

"""
Phenotypes from bin. pheno file
"""
def get_phenos(pheno_file):
    phenos = []
    with open(pheno_file, "r") as pheno_tab:
        line = pheno_tab.readline()
        line = line.rstrip('\n')
        line = line.split('\t')
        phenos = line
    return phenos

"""
Features from feature file
"""
def get_features(feature_file):
    features = []
    header = True
    with open(feature_file, "r") as feature_tab:
        for line in feature_tab:
            if header:
                header = False
                continue
            line = line.rstrip('\n')
            line = line.split('\t')
            features.append(line[0])
    return features

"""
Get samples from X or Y file, in_header = samples are in header 1st line), otherwise first column after header
"""
def get_samples(i_file, in_header=True):
    samples = []
    if in_header:
        with open(i_file, "r") as i_tab:
            line = i_tab.readline()
            line = line.rstrip('\n')
            line = line.split('\t')
            samples = line
    else:
        header = True
        with open(i_file, "r") as i_tab:
            for line in i_tab:
                if header:
                    header = False
                    continue
                line = line.rstrip('\n')
                line = line.split('\t')
                samples.append(line[0])
    return samples

"""
Pheno call stats
"""
def pheno_class_stats(pheno_file, pheno_name, samples=None):
    allowed_values = ['NA','1','0']
    pheno_class_count = {}
    header = True
    pheno_col = None
    with open(pheno_file, "r") as pheno:
        for line in pheno:
            line = line.rstrip('\n'); line = line.split('\t')
            if header: # header = phenotype names
                header = False
                pheno_col = [ i for i in range(0,len(line)) if line[i] == pheno_name ][0] + 1
                continue
            assert pheno_col is not None and pheno_col >= 1
            sampleID= line[0]
            pheno_v = str(line[pheno_col])
            assert pheno_v in allowed_values, "File %s, phenotype %s: got value %s, allowed are %s" %(pheno_file,pheno_name,pheno_v,', '.join(allowed_values))
            
            if sampleID not in samples:
                continue
            if pheno_v not in pheno_class_count:
                pheno_class_count[pheno_v] = 1
            else:
                pheno_class_count[pheno_v] += 1
    return pheno_class_count

"""
Pheno as list: values and sample IDs
"""
def pheno_as_list(pheno_file, pheno_name, ignore_miss=True, samples=None):
    allowed_values = ['NA','1','0']
    header = True
    pheno_col = None
    pheno_samples= []
    pheno_values = []
    with open(pheno_file, "r") as pheno:
        for line in pheno:
            line = line.rstrip('\n'); line = line.split('\t')
            if header: # header = phenotype names
                header = False
                pheno_col = [ i for i in range(0,len(line)) if line[i] == pheno_name ][0] + 1
                continue
            assert pheno_col is not None and pheno_col >= 1
            sampleID= line[0]
            pheno_v = str(line[pheno_col])
            assert pheno_v in allowed_values, "File %s, phenotype %s: got value %s, allowed are %s" %(pheno_file,pheno_name,pheno_v,', '.join(allowed_values))
            if sampleID not in samples:
                continue
            if ignore_miss and pheno_v == 'NA':
                continue
            pheno_samples.append(sampleID)
            pheno_values.append(pheno_v)
    assert len(pheno_values) == len(pheno_samples)
    return pheno_values, pheno_samples

"""
Sample intersection between X and Y with non-miss. value in Y and w.r.t. to optional include/exclude lists
"""
#def samples_intersect(pheno_file, geno_file, pheno_name, miss_char='NA', in_samples=None, ex_samples=None):
    #missing = 0
    #pheno_samples = set()
    #geno_samples  = set()
    ## geno: header = samples
    #with open(geno_file, "r") as geno:
        #line = geno.readline(); line = line.rstrip('\n'); line = line.split('\t')
        #geno_samples = set(line)
    ## pheno
    #header = True
    #pheno_col = None
    #with open(pheno_file, "r") as pheno:
        #for line in pheno:
            #line = line.rstrip('\n'); line = line.split('\t')
            #if header: # header = phenotype names
                #header = False
                #pheno_col = [ i for i in range(0,len(line)) if line[i] == pheno_name ][0] + 1
                #continue
            #assert pheno_col is not None and pheno_col >= 1
            #if line[0] in in_samples: sys.stdout.write('')
            #if line[pheno_col] != miss_char: # if pheno is not missing
                #pheno_samples.add(line[0])
            #else: missing += 1
    ## output: intersection
    #samples_i = pheno_samples.intersection(geno_samples)
    #if in_samples is not None and len(in_samples) > 0:
        #samples_i = in_samples.intersection(samples_i)
    #if ex_samples is not None:
        #samples_i = samples_i.difference(ex_samples)
    #sys.stdout.write("Column: %d\n" % pheno_col)
    #sys.stdout.write("Samples in Y (non-miss): %d\n" % len(pheno_samples))
    #sys.stdout.write("Missing in Y: %d\n" % missing)
    #sys.stdout.write("Intersection of non-miss in Y and X w.r.t. in/ex lists: %d\n" % len(samples_i))
    #return samples_i

# NOTE EIGENSTRAT specific stuff
# NOTE Read feature clusters: ID, cluster; no header; tab. sep.
def read_feature_cl(feature_cl_file):
    feature_cl = {}
    with open(feature_cl_file, "r") as feature_cls:
        for line in feature_cls:
            line = line.rstrip('\n')
            if line:
                line = line.split('\t')
                assert line[0] not in feature_cl
                feature_cl[line[0]] = line[1]
    return feature_cl

## NOTE sample intersection of the bin. pheno and geno file for ith phenotype
#def samples_intersect(pheno_file, geno_file, pheno_name, miss_char='NA', samples=None):
    #pheno_samples = set()
    #geno_samples = set()
    ## geno: header = samples
    #with open(geno_file, "r") as geno:
        #line = geno.readline(); line = line.rstrip('\n'); line = line.split('\t')
        #geno_samples = set(line)
    ## pheno
    #header = True
    #pheno_col = None
    #with open(pheno_file, "r") as pheno:
        #for line in pheno:
            #line = line.rstrip('\n'); line = line.split('\t')
            #if header: # header = phenotype names
                #header = False
                #pheno_col = [ i for i in range(0,len(line)) if line[i] == pheno_name ][0] + 1
                #continue
            #assert pheno_col is not None and pheno_col >= 1
            #if line[pheno_col] != miss_char: # if pheno is not missing
                #pheno_samples.add(line[0])
    ## output: intersection
    #if samples is not None:
        #return samples.intersection(pheno_samples.intersection(geno_samples))
    #else:
        #return pheno_samples.intersection(geno_samples)

# NOTE sample intersection of the bin. pheno and geno file for ith phenotype
# TODO obsolete -> mod. code in call_only_eigenstrat
def pheno_class_ratio(pheno_file, pheno_name, miss_char='NA', samples=None):
    pheno_classes = {}
    
    header = True
    pheno_col = None
    with open(pheno_file, "r") as pheno:
        for line in pheno:
            line = line.rstrip('\n'); line = line.split('\t')
            if header: # header = phenotype names
                header = False
                pheno_col = [ i for i in range(0,len(line)) if line[i] == pheno_name ][0] + 1
                continue
            assert pheno_col is not None and pheno_col >= 1
            sampleID= line[0]
            pheno_v = str(line[pheno_col])
            if sampleID not in samples:
                continue
            if pheno_v not in pheno_classes:
                pheno_classes[pheno_v] = 1
            else:
                pheno_classes[pheno_v] += 1
    total = sum(pheno_classes.values())
    max_ratio = 100.0 * float(max(pheno_classes.values())) / float(total)
    min_ratio = 100.0 * float(min(pheno_classes.values())) / float(total)
    
    pheno_classes['total'] = total
    pheno_classes['max_ratio'] = max_ratio
    pheno_classes['min_ratio'] = min_ratio
    return pheno_classes

# TODO obsolete -> mod. code in call_only_eigenstrat
def get_phenotypes(pheno_file):
    p_types = {}
    with open(pheno_file, "r") as pheno:
        line = pheno.readline()
        line = line.rstrip('\n')
        line = line.split('\t')
        p_types = dict.fromkeys(line,None)
    return p_types