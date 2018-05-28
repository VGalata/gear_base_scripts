#!/usr/bin/python

"""
TODO
    Based on nucleotide sequences
"""

"""
NOTE Resfams: http://www.nature.com/ismej/journal/v9/n1/full/ismej2014106a.html
Materials and Methods: Building of profile hidden Markov models (HMMs):
    Resfams AR family profile HMMs were built by
    (1) generating a multiple sequence alignment for each AR family (see Supplementary Methods)
    using MUSCLE (Edgar, 2004) v3.8.31 with default parameters and
    (2) training profile HMMs using the hmmbuild function of the HMMER3 (Finn et al., 2011) software package using default parameters.


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
Fastest possible (amino acids): -maxiters 1 -diags -sv -distance1 kbit20_3
Fastest possible (nucleotides): -maxiters 1 -diags

NOTE CMD {exe} -quiet -maxiters 2 -diags1 -in {ifile} -out {ofile}:
- chosen after calls with default parameters resulted in segmentation fault (A. baumannii, S. aureus)
- also some calls needed 2-4 hours
    - "-maxiters 2": http://www.drive5.com/muscle/manual/compromise.html
        "Iterations beyond 2 attempt refinement, which often results in a small improvement, at most."
    - "-diags1": http://www.drive5.com/muscle/manual/diagonals.html:
        "... speeds up the algorithm at the expense of some reduction in accuracy ..."
        "... reasonable strategy to enable diagonals in the first iteration but not the second ..."
- even with parameters above some nucl. alignments may need > 10h
    - http://www.drive5.com/muscle/manual/fastest.html
        for aminoacid seq.s: -maxiters 1 -diags -sv -distance1 kbit20_3
        fur nucl seq.s:      -maxiters 1 -diags


NOTE hmmbuild HMMER 3.1b2 (February 2015); http://hmmer.org/
Usage: hmmbuild [-options] <hmmfile_out> <msafile>

Basic options:
  -h     : show brief help on version and usage
  -n <s> : name the HMM <s>
  -o <f> : direct summary output to file <f>, not stdout
  -O <f> : resave annotated, possibly modified MSA to file <f>

Options for selecting alphabet rather than guessing it:
  --amino : input alignment is protein sequence data
  --dna   : input alignment is DNA sequence data
  --rna   : input alignment is RNA sequence data

Alternative model construction strategies:
  --fast           : assign cols w/ >= symfrac residues as consensus  [default]
  --hand           : manual construction (requires reference annotation)
  --symfrac <x>    : sets sym fraction controlling --fast construction  [0.5]
  --fragthresh <x> : if L <= x*alen, tag sequence as a fragment  [0.5]

Alternative relative sequence weighting strategies:
  --wpb     : Henikoff position-based weights  [default]
  --wgsc    : Gerstein/Sonnhammer/Chothia tree weights
  --wblosum : Henikoff simple filter weights
  --wnone   : don't do any relative weighting; set all to 1
  --wgiven  : use weights as given in MSA file
  --wid <x> : for --wblosum: set identity cutoff  [0.62]  (0<=x<=1)

Alternative effective sequence weighting strategies:
  --eent       : adjust eff seq # to achieve relative entropy target  [default]
  --eclust     : eff seq # is # of single linkage clusters
  --enone      : no effective seq # weighting: just use nseq
  --eset <x>   : set eff seq # for all models to <x>
  --ere <x>    : for --eent: set minimum rel entropy/position to <x>
  --esigma <x> : for --eent: set sigma param to <x>  [45.0]
  --eid <x>    : for --eclust: set fractional identity cutoff to <x>  [0.62]

Alternative prior strategies:
  --pnone       : don't use any prior; parameters are frequencies
  --plaplace    : use a Laplace +1 prior
  --popen <x>   : force gap open prob. (w/ --singlemx, aa default 0.02, nt 0.031)
  --pextend <x> : force gap extend prob. (w/ --singlemx, aa default 0.4, nt 0.75)

Handling single sequence inputs:
  --singlemx   : use substitution score matrix for single-sequence inputs
  --mx <s>     : substitution score matrix (built-in matrices, with --singlemx)
  --mxfile <f> : read substitution score matrix from file <f> (with --singlemx)

Control of E-value calibration:
  --EmL <n> : length of sequences for MSV Gumbel mu fit  [200]  (n>0)
  --EmN <n> : number of sequences for MSV Gumbel mu fit  [200]  (n>0)
  --EvL <n> : length of sequences for Viterbi Gumbel mu fit  [200]  (n>0)
  --EvN <n> : number of sequences for Viterbi Gumbel mu fit  [200]  (n>0)
  --EfL <n> : length of sequences for Forward exp tail tau fit  [100]  (n>0)
  --EfN <n> : number of sequences for Forward exp tail tau fit  [200]  (n>0)
  --Eft <x> : tail mass for Forward exponential tail tau fit  [0.04]  (0<x<1)

Other options:
  --cpu <n>          : number of parallel CPU workers for multithreads
  --stall            : arrest after start: for attaching debugger to process
  --informat <s>     : assert input alifile is in format <s> (no autodetect)
  --seed <n>         : set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
  --w_beta <x>       : tail mass at which window length is determined
  --w_length <n>     : window length 
  --maxinsertlen <n> : pretend all inserts are length <= <n>

NOTE HMMER and legal characters: http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf -> "Legal characters"
- Nucl: {ACGTU} + {NRY}
- Prot: 20 AA + {BXZ}

"""

__author__ = 'vgalata'

# NOTE IMPORTS
import sys
import os
import argparse
import time
import re
import glob
import shutil
import shlex, subprocess
from tqdm import tqdm # progress bar
import itertools
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO # parsing fasta files
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

__MUSCLE__='/home/vgalata/Programms/Muscle/muscle3.8.31_i86linux64'
__MUSCLE_PARAMS_PROT__='-maxiters 1 -diags -sv -distance1 kbit20_3' # http://www.drive5.com/muscle/manual/fastest.html
__MUSCLE_PARAMS_NUCL__='-maxiters 1 -diags' # http://www.drive5.com/muscle/manual/fastest.html
__MUSCLE_CMD__ = '{exe} -quiet {params} -in {ifile} -out {ofile}'
__HMMBUILD__='hmmbuild'
__HMMBUILD_CMD__='{exe} --seed 42 --cpu {cores} {alphabet} -n {name} -o {sfile} {ofile} {ifile}'
__HMM_LEGAL_NUCL_CHARS__=set(['A','C','G','T','U','N','R','Y'])
__HMM_LEGAL_PROT_CHARS__=set(['B','X','Z']).union( set(IUPAC.protein.letters) )
__TAR__='tar'
__LBZIP2_CMD__='{exe} --directory={sdir}{params} --format=pax --use=lbzip2 -c -f {obname}.tar.bz2 .'
__CAT_HMM_CMD__='cd {sdir}; cat {pattern} > {ofile}'


"""
"""

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

"""
"""

def timestamp():
    return time.strftime("%Y.%m.%d %H:%M:%S")

def scan_faa(fasta_file, seq_ids):
    if len(seq_ids) == 0:
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

"""
"""

def run_cmd(cmd, dry_run=False):
    if dry_run:
        return (cmd, '', 0)
    import shlex, subprocess
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    p_stdout = p.stdout.read().decode()
    p_comm   = p.communicate()[0]
    p_status = p.returncode
    return (cmd, p_stdout, p_status)

def run_muscle(ifile, params='', dry_run=False):
    cmd = __MUSCLE_CMD__.format(exe=__MUSCLE__, ifile=ifile, ofile="%s.msa" % ifile, params=params)
    return run_cmd(cmd, dry_run)

def run_hmmbuild(ifile, cores=1, alphabet='--amino', dry_run=False):
    bn1 = os.path.basename(ifile).replace('.msa','')
    bn2 = ifile.replace('.msa','')
    cmd = __HMMBUILD_CMD__.format(exe=__HMMBUILD__, ifile=ifile, alphabet=alphabet, cores=cores, name='centroid_%s' % bn1, sfile='%s.hmm_sum' % bn2, ofile='%s.hmm' % bn2)
    return run_cmd(cmd, dry_run)

def check_seq(s_id,s,alphabet):
    s_chars = set(str(s))
    if alphabet == 'prot':
        s_non_legal_chars = s_chars.difference(__HMM_LEGAL_PROT_CHARS__)
    elif alphabet == 'nucl':
        s_non_legal_chars = s_chars.difference(__HMM_LEGAL_NUCL_CHARS__)
    assert len( s_non_legal_chars ) == 0, 'SKIPPED: %s seq. of %s because contains %s' % (alphabet,s_id, ', '.join(list( s_non_legal_chars )))

"""
"""

def do_work(c_id__c_genes, odir, faa_files, ffn_files, cores=1): # NOTE had no fnn_files until 2016.09.28
    centroid = c_id__c_genes[0] # key
    c_genes  = c_id__c_genes[1] # values
    genes_faa = os.path.join(odir,'%s.faa' % centroid)
    #genes_ffn = os.path.join(odir,'%s.ffn' % centroid)
    genes_faa_msa="%s.msa" % genes_faa
    #genes_ffn_msa="%s.msa" % genes_ffn
    genes_faa_hmm="%s.hmm" % genes_faa
    #genes_ffn_hmm="%s.hmm" % genes_ffn
    
    if os.path.isfile(genes_faa_hmm) and os.path.isfile(genes_ffn_hmm):
        sys.stdout.write('SKIP: %s: HMMs exist: %s, %s\n' % (centroid, genes_faa_hmm, genes_ffn_hmm))
        return
    
    # protein sequences
    prot_seqs = {}
    # FAA (protein sequences)
    pool        = Pool(cores) # create pool for parallel computing
    pool_iter   = itertools.product(faa_files, [c_genes]) # iterable for the pool
    results     = pool.starmap( scan_faa , pool_iter ) # perform parallel computing
    pool.close(); pool.join() # wait until all finished and close the pool
    for records in results:
        assert all( [ r_id not in prot_seqs for r_id in records.keys() ] ) # TEST
        prot_seqs.update( records )
    # FFN (translate nucl. sequences)
    remaining = c_genes.difference( prot_seqs.keys() )
    if len(remaining) > 0:
        sys.stdout.write('\tINFO: for %d genes need to scan ffn files\n' % len(remaining))
        pool            = Pool(cores) # create pool for parallel computing
        pool_iter       = itertools.product(ffn_files, [remaining], [True]) # iterable for the pool NOTE added True for translation 2016.09.28
        results         = pool.starmap( scan_ffn , pool_iter ) # perform parallel computing NOTE had scan_faa until 2016.09.28
        pool.close(); pool.join() # wait until all finished and close the pool
        for records in results:
            assert all( [ r_id not in prot_seqs for r_id in records.keys() ] ) # TEST
            prot_seqs.update( records )
    assert all([ r_id in c_genes for r_id in prot_seqs.keys() ]), ';'.join( list(set(prot_seqs.keys()).difference( c_genes )) ) # TEST NOTE changed to prot_seqs - c_genes since 2016.09.28
    assert all([ r_id in prot_seqs.keys() for r_id in c_genes ]), ';'.join( list(c_genes.difference( prot_seqs.keys() )) ) # TEST NOTE New, since 2016.09.28
    
    ## nucl. sequences
    #nucl_seqs = {}
    #pool        = Pool(cores) # create pool for parallel computing
    #pool_iter   = itertools.product(ffn_files, [c_genes], [False]) # iterable for the pool
    #results     = pool.starmap( scan_ffn , pool_iter ) # perform parallel computing
    #pool.close(); pool.join() # wait until all finished and close the pool
    #for records in results:
        #assert all( [ r_id not in nucl_seqs for r_id in records.keys() ] ) # TEST
        #nucl_seqs.update( records )
    #assert all( [r_id in c_genes for r_id in nucl_seqs.keys()] ), ';'.join( list(c_genes.difference( nucl_seqs.keys() )) ) # TEST
    
    # write FASTA
    with open(genes_faa, "w") as f:
        for c_gene in c_genes:
            check_seq(c_gene,prot_seqs[ c_gene ].seq,'prot')
            SeqIO.write(prot_seqs[ c_gene ],f, "fasta")
    #with open(genes_ffn, "w") as f:
        #for c_gene in c_genes:
            #check_seq(c_gene,nucl_seqs[ c_gene ].seq,'nucl')
            #SeqIO.write(nucl_seqs[ c_gene ],f, "fasta")
    
    # run MUSCLE
    cmd, cmd_stdout, cmd_status = run_muscle(ifile=genes_faa, params=__MUSCLE_PARAMS_PROT__)
    assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout) # TEST
    #cmd, cmd_stdout, cmd_status = run_muscle(ifile=genes_ffn)
    #assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout) # TEST
    # build HMM
    cmd, cmd_stdout, cmd_status = run_hmmbuild(ifile=genes_faa_msa, alphabet='--amino')
    assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout) # TEST
    #cmd, cmd_stdout, cmd_status = run_hmmbuild(ifile=genes_ffn_msa, alphabet='--dna')
    #assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout) # TEST
    # rm fasta
    os.remove(genes_faa)
    #os.remove(genes_ffn)

def compress_MSA(sdir, obname, dry_run=False):
    cmd = __LBZIP2_CMD__.format(exe=__TAR__, sdir=sdir, obname=obname, params=' --exclude=\'*.hmm*\'')
    return run_cmd(cmd, dry_run)

def collect_hmm(sdir, ofile, pattern, dry_run=False):
    cmd = __CAT_HMM_CMD__.format(sdir=sdir, ofile=ofile, pattern=pattern)
    return run_cmd(cmd, dry_run)

# NOTE MAIN
if __name__ == "__main__":
    # NOTE Args
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    # required
    parser.add_argument('--ref',        '-f', help='reference FASTA file created by Roary', required=True)
    parser.add_argument('--pam',        '-p', help='presence-absence matrix created by Roary', required=True)
    parser.add_argument('--prokka',           help='', required=True)
    parser.add_argument('--odir',       '-o', help='output file', required=True)
    # optional
    parser.add_argument('--ncol',       '-n', help='in prokka files table nth column contains the file names (default 1)', default=1, type=int)
    parser.add_argument('--cores',      '-c', help='cores to use', default=1, type=int)
    parser.add_argument('--job_cores',  '-j', help='cores per job', default=1, type=int)
    parser.add_argument('--test',       '-t', help="Test mode", action="store_true")
    parser.add_argument('--verbose',    '-v', help="Verbosity", action="store_true")
    args = parser.parse_args()
    
    if args.verbose:
        sys.stdout.write('MUSCLE: %s: %s with\n\tprot: %s\n\tnucl: %s\n' % (__MUSCLE__,__MUSCLE_CMD__,__MUSCLE_PARAMS_PROT__,__MUSCLE_PARAMS_NUCL__))
        sys.stdout.write('HMMER : %s: %s\n' % (__HMMBUILD__,__HMMBUILD_CMD__))
        sys.stdout.write("%s: Start.\n" % timestamp())
    
    # NOTE Init
    test_count = 3
    if not os.path.exists(os.path.join(args.odir,'tmp')):
        os.makedirs(os.path.join(args.odir,'tmp'))
    
    # Prokka files
    faa_files = []
    ffn_files = []
    with open(args.prokka, 'r') as prokka_files:
        for line in prokka_files:
            if line:
                f = line.rstrip('\n').split('\t')[args.ncol-1]
                faa_files.append( "%s.faa" % os.path.splitext(f)[0] )
                ffn_files.append( "%s.ffn" % os.path.splitext(f)[0] )
    
    # FASTA
    seqs = SeqIO.parse(open(args.ref),'fasta')
    centroid_dict = {}
    for centroid in [ fasta_entry.description for fasta_entry in seqs ]:
        centroidID = centroid.split(' ')[0]
        centroidName = ' '.join(centroid.split(' ')[1:]) # Name given by Roary
        assert centroidName not in centroid_dict, "Centroid %s already in dict with ID=%s" % (centroidName,centroid_dict[centroidName])
        centroid_dict[centroidName] = centroidID
        if args.test and len(centroid_dict) == test_count: # TEST
            break
    if args.verbose:
        sys.stdout.write("%s: Centroids: %d\n" % (timestamp(),len(centroid_dict)))
    
    # Presence/absence matrix -> get all clustered genes
    centroids = dict.fromkeys(centroid_dict.values(),None)
    all_genes = set()
    with open(args.pam,'r') as ipam:
        header = True
        for line in ipam:
            line = line.rstrip('\r\n') # remove new line
            line = line.split('","') # tab. is "," sep. but some entries may also have ",", additionally each entry is quoted
            line = [re.sub(r'"', '', entry) for entry in line] # remove remaining quotes
            
            if header:
                header = False
                continue
            
            centroidName    = line[0] # centroid name/ID assigned by Roary
            clustered_genes = line[14:] # get clustered genes
            
            if args.test and centroidName not in centroid_dict: # TEST
                continue
            
            assert centroidName in centroid_dict, '%s not found' % centroidName
            centroidID = centroid_dict[centroidName]
            
            centroids[centroidID] = set()
            for cl_genes in clustered_genes:
                if cl_genes == "": # this sample has no gene belonging to that cluster
                    continue
                for cl_gene in cl_genes.split('\t'):
                    assert cl_gene not in centroids[centroidID], '%s already found' % cl_gene
                    assert cl_gene not in all_genes, '%s already found' % cl_gene
                    centroids[centroidID].add( cl_gene )
                    all_genes.add( cl_gene )
    if args.verbose:
        sys.stdout.write("%s: Genes: %d\n" % (timestamp(), len(all_genes)))
    
    # NOTE do work: extract protein sequences, create MSA, build HMM
    pool            = MyPool(int(args.cores/args.job_cores)) # create pool for parallel computing
    pool_iter       = itertools.product(centroids.items(), ['%s/tmp' % args.odir], [faa_files], [ffn_files], [args.job_cores]) # iterable for the pool NOTE had no ffn_files until 2016.09.28
    results         = pool.starmap( do_work , pool_iter ) # perform parallel computing
    pool.close(); pool.join() # wait until all finished and close the pool
    
    if args.verbose:
        sys.stdout.write("%s: Created HMM models\n" % timestamp() )
    
    # NOTE Clean up
    # archive *.msa
    if args.verbose:
        sys.stdout.write("%s: Creating MSA archive\n" % timestamp() )
    cmd, cmd_stdout, cmd_status = compress_MSA(sdir='%s/tmp' % args.odir, obname=os.path.join(args.odir,'MSA'))
    if args.verbose:
        sys.stdout.write('\t%s\n' % cmd)
    assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout) # TEST
    # sync *.hmm
    if args.verbose:
        sys.stdout.write("%s: Moving *.hmm files\n" % timestamp() )
    cmd, cmd_stdout, cmd_status = collect_hmm(sdir='%s/tmp/' % args.odir, ofile=os.path.join(args.odir,'pangenome.faa.hmm'), pattern='*.faa.hmm')
    if args.verbose:
        sys.stdout.write('\t%s\n' % cmd)
    assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout) # TEST
    #cmd, cmd_stdout, cmd_status = collect_hmm(sdir='%s/tmp/' % args.odir, ofile=os.path.join(args.odir,'pangenome.ffn.hmm'), pattern='*.ffn.hmm')
    #if args.verbose:
        #sys.stdout.write('\t%s\n' % cmd)
    #assert cmd_status == 0, '%s: %d: %s' % (cmd, cmd_status, cmd_stdout) # TEST
    # rm tmp dir.
    if args.verbose:
        sys.stdout.write("%s: Removing tmp dir.\n" % timestamp() )
    shutil.rmtree(os.path.join(args.odir, 'tmp'))
    
    if args.verbose:
        sys.stdout.write("%s: Finished.\n" % timestamp())
