#!/usr/bin/python

"""
NOTE MUSCLE v3.8.31 by Robert C. Edgar
Basic usage
muscle -in <inputfile> -out <outputfile>
Common options (for a complete list please see the User Guide):
    -in <inputfile>    Input file in FASTA format (default stdin)
    -out <outputfile>  Output alignment in FASTA format (default stdout)
    -diags             Find diagonals (faster for similar sequences)
    -maxiters <n>      Maximum number of iterations (integer, default 16)
    -maxhours <h>      Maximum time to iterate in hours (default no limit)
    -html              Write output in HTML format (default FASTA)
    -msf               Write output in GCG MSF format (default FASTA)
    -clw               Write output in CLUSTALW format (default FASTA)
    -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
    -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
    -quiet             Do not write progress messages to stderr
    -version           Display version information and exit
Without refinement (very fast, avg accuracy similar to T-Coffee): -maxiters 2
"""

__author__ = 'vgalata'

# NOTE IMPORTS
import sys
import os
import argparse
import re
import glob
import shutil
import shlex, subprocess
from tqdm import tqdm # progress bar
import itertools
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool
import logging
import pandas

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO # parsing fasta files
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

__MUSCLE__='/home/vgalata/Programms/Muscle/muscle3.8.31_i86linux64'
__MUSCLE_PARAMS_PROT__='-maxiters 1 -diags -sv -distance1 kbit20_3' # http://www.drive5.com/muscle/manual/fastest.html
__MUSCLE_CMD__ = '{exe} -quiet {params} -in {ifile} -out {ofile}'

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

# sample ID from gene ID
def sample_id_from_gene_id(g_id):
    # pattern: <SampleID>_<geneID>
    return g_id.split('_')[0]

# sample ID from GFF file name
def sample_id_from_prokka_file(f):
    # pattern: <path>/<SampleID>.gff
    # ID is the basename w/o extension
    return os.path.splitext(os.path.basename(f))[0]

def get_msa_file_name(odir, k):
    return os.path.join(odir, 'msa_%s' % str(k))

def scan_faa(fasta_file, seq_ids):
    if len(seq_ids) == 0:
        return {}
    if sample_id_from_prokka_file(fasta_file) not in [sample_id_from_gene_id(g_id) for g_id in seq_ids]:
        return {}
    records = {}
    seqs = SeqIO.parse(open(fasta_file),'fasta')
    for entry in seqs:
        if entry.id in seq_ids:
            assert entry.id not in records, "%s already added" % entry.id
            entry.description = ''
            records[entry.id] = entry
    return records

def scan_ffn(fasta_file, seq_ids, translate_seq=True):
    if len(seq_ids) == 0:
        return {}
    if sample_id_from_prokka_file(fasta_file) not in [sample_id_from_gene_id(g_id) for g_id in seq_ids]:
        return {}
    records = {}
    seqs = SeqIO.parse(open(fasta_file),'fasta')
    for entry in seqs:
        if entry.id in seq_ids:
            if translate_seq:
                # translate
                entry_prot = entry.seq.translate(table=11, stop_symbol='*',to_stop=False, cds=False)
                # remove stop codon
                if entry_prot[len(entry_prot)-1] == '*':
                    entry_prot = entry_prot[:(len(entry_prot)-1)]
                entry_prot = str(entry_prot)
                # create record
                assert entry.id not in records, "%s already added" % entry.id
                records[entry.id] = SeqRecord(Seq(entry_prot, IUPAC.protein), id=centroid, description='')
            else:
                entry.description = ''
                records[entry.id] = entry
    return records

def run_cmd(cmd, dry_run=False):
    if dry_run:
        return (cmd, '', 0)
    import shlex, subprocess
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    p_stdout = p.stdout.read().decode()
    p_comm   = p.communicate()[0]
    p_status = p.returncode
    return (cmd, p_stdout, p_status)

def run_muscle(ifile, ofile, params='', dry_run=False):
    cmd = __MUSCLE_CMD__.format(exe=__MUSCLE__, ifile=ifile, ofile=ofile, params=params)
    return run_cmd(cmd, dry_run)

def do_work(egene, odir, faa_files, ffn_files, cores=1):
    # get best hits for given essential gene
    set_genes = set(list(hits.loc[hits['query_name']==egene]['target_name']))
    logger.info('Essential gene %s: %d genes' % (egene, len(set_genes)))

    # output files: FAA with prot. sequences and MSA
    genes_faa = os.path.join(odir,'gene_set_%s.faa' % egene)
    genes_faa_msa = "%s.msa" % os.path.splitext(genes_faa)[0]

    # protein sequences
    prot_seqs = {}

    # FAA (protein sequences)
    pool        = Pool(cores) # create pool for parallel computing
    pool_iter   = itertools.product(faa_files, [set_genes])
    results     = pool.starmap( scan_faa , pool_iter )
    pool.close(); pool.join()
    for records in results:
        assert all( [ r_id not in prot_seqs for r_id in records.keys() ] )
        prot_seqs.update(records)

    # FFN (translate nucl. sequences)
    remaining = set_genes.difference( prot_seqs.keys() )
    if len(remaining) > 0:
        logger.info('FFNs required for %d genes for SET %s' % (len(remaining), egene))
        pool            = Pool(cores) # create pool for parallel computing
        pool_iter       = itertools.product(ffn_files, [remaining], [True])
        results         = pool.starmap( scan_ffn , pool_iter )
        pool.close(); pool.join()
        for records in results:
            assert all( [ r_id not in prot_seqs for r_id in records.keys() ] ) # TEST
            prot_seqs.update( records )

    # Check if all found
    assert all([ r_id in set_genes for r_id in prot_seqs.keys() ]), ';'.join( list(set(prot_seqs.keys()).difference( set_genes )) )
    assert all([ r_id in prot_seqs.keys() for r_id in set_genes ]), ';'.join( list(set_genes.difference( prot_seqs.keys() )) )

    # write FASTA
    with open(genes_faa, "w") as ofile:
        for gene in set_genes:
            SeqIO.write(prot_seqs[gene], ofile, "fasta")

    # run MUSCLE
    cmd, cmd_stdout, cmd_status = run_muscle(ifile=genes_faa, ofile=genes_faa_msa, params=__MUSCLE_PARAMS_PROT__)
    assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout)
    return

# NOTE MAIN
if __name__ == "__main__":
    # Logger
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Args
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

    parser.add_argument('--hits', help='Hits: samples x genes', required=True)
    parser.add_argument('--egenes', help='Seleced genes (list)', required=True)
    parser.add_argument('--samples', help='Seleced samples (list)', required=True)
    parser.add_argument('--prokka', help='Tables with PROKKA files used to build pan-genomes', required=True)
    parser.add_argument('--odir', help='Output directory', required=True)

    parser.add_argument('--cores',      '-c', help='cores to use', default=1, type=int)
    parser.add_argument('--job_cores',  '-j', help='cores per job', default=1, type=int)
    parser.add_argument('--verbose',    '-v', help="Verbosity", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logger.info('MUSCLE: %s: %s with\n\tprot: %s' % (__MUSCLE__,__MUSCLE_CMD__,__MUSCLE_PARAMS_PROT__))

    # Prokka files
    faa_files = []
    ffn_files = []
    with open(args.prokka, 'r') as ifile:
        for line in ifile:
            f = line.rstrip('\n')
            f_b = os.path.basename(f)
            f_d = os.path.dirname(f)
            faa_files.append( "%s.faa" % os.path.splitext(f)[0] )
            ffn_files.append( "%s.ffn" % os.path.splitext(f)[0] )

    # Data: hits, selected essential genes and samples
    hits = pandas.read_csv(filepath_or_buffer=args.hits, sep='\t', header=0, index_col=None)
    logger.info('Number of best hits: %d' % hits.shape[0])

    egenes = list(pandas.read_csv(filepath_or_buffer=args.egenes, sep='\t', header=None, index_col=None)[0])
    logger.info('Number of selected essential genes: %d' % len(egenes))

    samples = list(pandas.read_csv(filepath_or_buffer=args.samples, sep='\t', header=None, index_col=None)[0])
    logger.info('Number of selected samples: %d' % len(samples))

    # Reduce hits
    hits = hits.loc[(hits['sample'].isin(samples)) & (hits['query_name'].isin(egenes))]
    logger.info('Number of best hits after filtering: %d' % hits.shape[0])

    # Extract protein sequences, create MSA
    if args.verbose: logger.info('Extract AA seq.s and perform MSA...')
    pool            = MyPool(int(args.cores/args.job_cores)) # create pool for parallel computing
    pool_iter       = itertools.product(egenes, [args.odir], [faa_files], [ffn_files], [args.job_cores])
    results         = pool.starmap( do_work , pool_iter ) # perform parallel computing
    pool.close(); pool.join() # wait until all finished and close the pool

    if args.verbose: logger.info('Done.')
