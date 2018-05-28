#!/usr/bin/python

__author__ = 'vgalata'

# NOTE FORMATS
"""
GFF 3 (http://www.sequenceontology.org/gff3.shtml)
Undefined fields are replaced with the "." character, as described in the original GFF spec.
Column 1: "seqID"
  The ID of the landmark used to establish the coordinate system for the current feature. [...]
Column 2: "source"
  The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature. [...]
Column 3: "type"
  The type of the feature (previously called the "method").
Columns 4 & 5: "start" and "end"
  The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative to the landmark given in column one. [...]
Column 6: "score"
Column 7: "strand"
  The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded.
  In addition, ? can be used for features whose strandedness is relevant, but unknown.
Column 8: "phase"
  For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. [...]
  The phase is REQUIRED for all CDS features.
Column 9: "attributes"
  A list of feature attributes in the format tag=value [...]
  Separated by ";"

VCF
Header/comment lines start with "#"
Column names (last header/comment line): CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [Sample ID 1] ... [Sample ID n]
Columns are sep. by "\t"
Fields in "INFO" are separated by ";"
snpEff annotations are in column "INFO" starting with "ANN=", alleles are sep. by ",", fields of one allele are sep. by "|"
    Example: ANN=A|missense_variant|MODERATE|rimP|C3830-210_00002|transcript|C3830-210_00002|protein_coding|1/1|c.340G>A|p.Gly114Ser|340/525|340/525|114/174||
"""

# NOTE IMPORTS
import sys
import os
import argparse
import time
import re
import subprocess
import Bio
from Bio import SeqIO # parsing fasta files


# NOTE METHODS

"""
Get annotation from GFF file by feature ID
"""
def annot_from_gff(gff_file, ids, attrs=['gene','product']):
    annot = {}
    IDs = set(ids) # create set of IDs
    with open(gff_file, 'r') as gff:
        for line in gff:
            line = line.rstrip('\n') # remove new line
            if not line or re.match("^##",line): # empty, header or FASTA
                if re.match("^##FASTA",line): # FASTA -> stop
                    break
                else:
                    continue
            try:
                f_seq_ID, f_source, f_type, f_start, f_end, f_score, f_strand, f_phase, f_attrs = line.split('\t') # columns are separated by tabs
                f_attrs = f_attrs.split(';')
            except ValueError:
                sys.stdout.write("Value error in line: %s" % line)
            except AttributeError:
                sys.stdout.write("Attribute error in line: %s" % line)
            f_ID = f_attrs[0] # 1st attr is the ID as "ID=<id>"
            for ID in IDs: # check each ID
                if f_ID == ("ID=%s" % ID): # Found ID
                    assert ID not in annot, "ID %s already in annotation dict"
                    annot[ID] = dict.fromkeys(attrs,"NA") # empty dict: NA for each asked attribute
                    for ID_attr in f_attrs:
                        for attr in attrs:
                            if re.match("^%s=.*" % attr,ID_attr):
                                annot[ID][attr]    = ID_attr.replace("%s=" % attr,"")
                    break
    return annot

"""
Get annotations from VCF
"""
def get_annot_from_vcf(vcf_file, ids=None):
    attrs=['chrom','pos','ref','alt','gene','DNA_change','AA_change']
    if ids is not None:
        ids = set(ids)                      # be sure to have unique IDs
        assert len(ids) > 0
        annot = dict.fromkeys(ids, None)    # annotation dictionary
    else:
        annot = {}
    count = 1                               # current variant ID
    with open(vcf_file,'r') as vcf_i:
        for line in vcf_i:
            line = line.rstrip('\n')        # remove new line
            if not line: continue           # empty line -> skip
            if re.match('^#',line):         # header lines -> skip
                continue
            line_s = line.split('\t')       # split into fields
            info   = line_s[7].split(';')   # split info into fields
            info   = [i.replace('ANN=','') for i in info if re.match('^ANN',i)][0] # get field with annot
            info   = info.split(',')        # split into alleles
            if (ids is not None and count in ids) or ids is None:
                annot[count] = dict.fromkeys(attrs, None)
                annot[count]['chrom']   = line_s[0]
                annot[count]['pos']     = line_s[1]
                annot[count]['ref']     = line_s[3]
                annot[count]['alt']     = line_s[4]
                annot[count]['gene']    = info[0].split('|')[4]
                annot[count]['DNA_change']  = '|'.join([ s.split('|')[9]  for s in info ])
                annot[count]['AA_change']   = '|'.join([ s.split('|')[10] for s in info ])
            count += 1                  # increase counter
    if ids is not None:
        assert len(annot) == len(ids), "There are %d unique IDs but %d were found" % (len(ids),len(annot))
    return annot