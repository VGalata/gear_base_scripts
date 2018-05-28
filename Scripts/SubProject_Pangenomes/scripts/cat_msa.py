#!/usr/bin/python

__author__ = 'vgalata'

# NOTE IMPORTS
import sys
import os
import argparse
import re
from tqdm import tqdm
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# sample ID from gene ID
def sample_id_from_gene_id(g_id):
    # pattern: <SampleID>_<geneID>
    return g_id.split('_')[0]

# MAIN
if __name__ == "__main__":
    # Args
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--msa', help='MSA files', required=True, nargs='+')
    parser.add_argument('--samples', help='Sample IDs', required=True)
    parser.add_argument('--ofile', help='Output file', required=True)
    parser.add_argument('--verbose',    '-v', help="Verbosity", action="store_true")
    args = parser.parse_args()
    args.msa = sorted(args.msa)

    # Logger
    # logging.basicConfig(stream=sys.stdout, format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.INFO)

    logFormatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger()

    fileHandler = logging.FileHandler("%s.log" % args.ofile)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    # if args.verbose:
    #     consoleHandler = logging.StreamHandler()
    #     consoleHandler.setFormatter(logFormatter)
    #     logger.addHandler(consoleHandler)

    # Read in samples
    samples = set()
    with open(args.samples, 'r') as ifile:
        for line in ifile:
            samples.add(line.rstrip('\n').strip())
    if args.verbose:
        logger.info('Samples: %d' % len(samples))

    # Collect MSAs
    msa_cat = dict.fromkeys(samples, None)
    for sample_id in samples:
        msa_cat[sample_id] = SeqRecord(Seq('', IUPAC.protein), id=sample_id, description='')
    for msa_file in tqdm(args.msa):
        msa = SeqIO.parse(open(msa_file),'fasta')
        msa_samples = set()
        msa_len = 0
        for entry in tqdm(msa):
            # set seq. ID to sample ID, i.e. remove gene ID part
            entry.id = sample_id_from_gene_id(entry.id)
            # log
            logger.info('%s: %s: %d' % (msa_file, entry.id, len(entry.seq)))
            # set MSA length if first seq.
            if msa_len == 0:
                msa_len = len(entry.seq)
            # checks
            assert entry.id not in msa_samples
            assert msa_len == len(entry.seq)
            # concat
            msa_cat[entry.id].seq += entry.seq
            msa_samples.add(entry.id)
        # check: max. MSA length = min. MSA length + current length
        total_msa_len = [len(s.seq) for s in msa_cat.values()]
        assert min(total_msa_len) == max(total_msa_len) or min(total_msa_len) == (max(total_msa_len) - msa_len)
        # for samples without MSA add empty MSA
        for sample_id in tqdm(samples.difference(msa_samples)):
            # check
            assert len(msa_cat[sample_id]) == max(total_msa_len) - msa_len
            # log
            logger.info('%s: %s: empty seq.' % (msa_file, sample_id))
            # add empty MSA seq.
            msa_cat[sample_id].seq += Seq('-' * msa_len, IUPAC.protein)

    # Write to file
    with open(args.ofile, 'w') as ofile:
        for sample_id in sorted(samples):
            SeqIO.write(msa_cat[sample_id], ofile, "fasta")
    logger.info('Done.')
