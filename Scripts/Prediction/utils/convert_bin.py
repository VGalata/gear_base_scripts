#!/bin/python

# NOTE Imports
import sys
import os
import argparse
import time
import re
from collections import OrderedDict

from utils.my_utils import read_id_list

def bin_to_plink_ped_map(geno_file, map_ofile, ped_ofile, family="Fam", sample_file=None, feature_file=None, feature_cl_file=None, miss='NA', verbose=True):
    if verbose: sys.stdout.write("Bin. feature mat. to PLINK PED/MAP")
    samples = None
    if sample_file is not None:
        samples = set(read_id_list(sample_file))
        if verbose: sys.stdout.write("Sample list contains %d unique IDs\n" % len(samples))
    features = None
    if feature_file is not None:
        features = set(read_id_list(feature_file))
        if verbose: sys.stdout.write("Feature list contains %d unique IDs\n" % len(features))
    feature_cl = None
    if feature_cl_file is not None:
        feature_cl = read_feature_cl(feature_cl_file)
    # MAP
    header = True
    count  = 1
    with open(geno_file, "r") as geno, open(map_ofile, "w") as map_o:
        if verbose: sys.stdout.write("\tOutput: %s\n" % map_ofile)
        for line in geno:
            if header:
                header = False
                continue
            line   = line.rstrip('\n')
            snpID  = line.split('\t')[0]
            snpChr = "0"
            snpDis = "0"
            snpPos = str(count)
            if features is not None and snpID not in features:
                continue
            if feature_cl is not None:
                snpChr = feature_cl[snpID]
            map_o.write("\t".join([snpChr,snpID,snpDis,snpPos]) + "\n")
            #map_o.flush()
            count += 1
    # PED
    header = True
    sample_dict = None
    sample_IDs  = None
    with open(geno_file, "r") as geno:
        if verbose: sys.stdout.write("\tOutput: %s\n" % ped_ofile)
        for line in geno:
            line = line.rstrip('\n'); line = line.split('\t')
            if header:
                header = False
                sample_IDs  = line
                sample_dict = dict.fromkeys(sample_IDs, "")
                continue
            snpID = line[0]
            for i in range(1,len(line)):
                sampleID = sample_IDs[i-1]
                if samples is not None and sampleID not in samples:
                    if sampleID in sample_dict: del sample_dict[sampleID]
                    continue
                if features is not None and snpID not in features:
                    continue
                sampleAl = line[i]
                if sampleAl == "1":
                    sampleAl = "A A"
                elif sampleAl == "0":
                    sampleAl = "G G"
                elif sampleAl == miss:
                    sampleAl = "0 0"
                if sample_dict[sampleID] == "":
                    sample_dict[sampleID] = sampleAl
                else:
                    sample_dict[sampleID] += ("\t" + sampleAl)
    with  open(ped_ofile, "w") as ped_o:
        for sampleID, sampleAl in sample_dict.items():
            ped_o.write("\t".join([family,sampleID,"0","0","0","0",sampleAl]) + "\n")
            #ped_o.flush()

def bin_to_plink_pheno(pheno_file, pheno_name, pheno_ofile, family="Fam", sample_file=None, miss='NA', verbose=True):
    if verbose: sys.stdout.write("Bin. pheno mat. to PLINK alt. pheno")
    samples = None
    if sample_file is not None:
        samples = set(read_id_list(sample_file))
        if verbose: sys.stdout.write("Sample list contains %d unique IDs\n" % len(samples))
    # Pheno
    header = True
    pheno_col = 0
    with open(pheno_file, "r") as pheno, open(pheno_ofile, "w") as pheno_o:
        sys.stdout.write("\tOutput: %s in %s\n" % (pheno_name, pheno_ofile))
        for line in pheno:
            line = line.rstrip('\n'); line = line.split('\t')
            if header:
                header = False
                for i in range(0,len(line)):
                    if line[i] == pheno_name:
                        pheno_col = i+1 # +1 because 1s field is sampleID
                        break
                continue
            sampleID = line[0]
            if samples is not None and sampleID not in samples:
                continue
            samplePh = line[pheno_col]
            if samplePh == miss:
                samplePh = '-9'
            pheno_o.write("\t".join([family,sampleID,samplePh])+"\n")
            #pheno_o.flush()

"""
NOTE Binary matrices to EIGENSTRAT (smart) format
"""
def bin_to_eig_geno_snp(geno_file, geno_ofile, snp_ofile, sample_file=None, feature_file=None, miss='NA'):
    samples = features = None
    if sample_file is not None:
        samples = set(read_id_list(sample_file))
    if feature_file is not None:
        features = set(read_id_list(feature_file))
    # SNP
    header = True
    count  = 1
    found_features = 0
    with open(geno_file, "r") as geno, open(snp_ofile, "w") as snp_o:
        for line in geno:
            if header:
                header = False
                continue
            line = line.rstrip('\n')
            snpID = line.split('\t')[0]
            if features is not None and snpID not in features:
                continue
            found_features += 1
            snp_o.write("\t".join(["rs"+snpID,"1","0.0",str(count),"",""]) + "\n") # EIGENSTRAT expects "rs<snpID>"
            count += 1
    assert (features is None) or (found_features== len(features)), "%s: expected %d, found %d" % (snp_file,len(features),found_features)
    # Geno
    header = True
    sampleIDs = None
    found_features = 0
    with open(geno_file, "r") as geno, open(geno_ofile, "w") as geno_o:
        for line in geno:
            line = line.rstrip('\n')
            line = line.split('\t')
            if header:
                header = False
                sampleIDs = line
                continue
            snpID = line[0]
            if features is not None and snpID not in features: continue
            found_features += 1
            found_samples   = 0
            for i in range(1,len(line)):
                sampleID = sampleIDs[i-1]
                if (samples is not None and sampleID not in samples) or (features is not None and snpID not in features):
                    continue
                found_samples += 1
                sampleAl = line[i]
                if sampleAl == miss: continue #sampleAl = "9"
                else: geno_o.write("\t".join(["rs"+snpID,sampleID,sampleAl]) + "\n")
            assert (samples is None) or (found_samples == len(samples)), "%s: expected %d, found %d" % (geno_file,len(saples),found_samples)
    assert (features is None) or (found_features== len(features)), "%s: expected %d, found %d" % (geno_file,len(features),found_features)

def bin_to_eig_pheno(pheno_file, pheno_name, pheno_ofile, sample_file=None, miss='NA'):
    samples = None
    if sample_file is not None:
        samples = set(read_id_list(sample_file))
    # Pheno
    header = True
    pheno_col = 0
    with open(pheno_file, "r") as pheno, open(pheno_ofile, "w") as pheno_o:
        for line in pheno:
            line = line.rstrip('\n')
            line = line.split('\t')
            if header:
                header = False
                for i in range(0,len(line)):
                    if line[i] == pheno_name:
                        pheno_col = i+1 # +1 because 1st field is sampleID
                        break
                continue
            assert pheno_col > 0
            sampleID = line[0]
            if samples is not None and sampleID not in samples:
                continue
            samplePh = line[pheno_col]
            if samplePh == '1': samplePh = 'Case'
            elif samplePh == '0': samplePh = 'Control'
            elif samplePh == miss: samplePh = 'Ignore'
            else: sys.exit('Unknown phenotype: %s: %s: %s' % (pheno_name,sampleID,samplePh))
            pheno_o.write("\t".join([sampleID,'U',samplePh])+"\n") #; pheno_o.flush()