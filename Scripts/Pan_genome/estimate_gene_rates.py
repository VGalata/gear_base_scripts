#!/usr/bin/python

__author__ = 'vgalata'

# NOTE Imports
import sys
import os
import argparse
import time
import re
import numpy
import itertools
from multiprocessing import Pool

def estimate_numbers(perm_i,m,cfreq=90):
    # permute columns (samples/genomes)
    m = m[:,numpy.random.permutation(m.shape[1])]
    #print(m.astype(int))
    # init: no gene present
    #genes       = numpy.zeros((m.shape[0]), dtype=bool) # whether found at least once
    genes   = numpy.zeros((m.shape[0]), dtype=int) # how many times the genes was found
    num     = numpy.zeros((m.shape[1],3+len(cfreq)), dtype=int) # matrix for results
    # for each column (sample/genome)
    for i in range(m.shape[1]):
        # new = not found so far but in this genome
        num_new     = sum(numpy.logical_and( numpy.logical_not( genes.astype(bool) ) , m[:,i] ))
        # update
        genes       = genes + m[:,i]
        # all/core/unique
        num_all     = sum( genes.astype(bool) )
        #num_core    = sum( [ (100.0 * float(g_n) / float(i+1)) >= cfreq for g_n in genes ] ) # core = found in at least x% of i
        num_uniq    = sum( [ g_n == 1 for g_n in genes ] ) # unique = found only once
        # save
        #num[i,0:3] = numpy.array([num_all,num_new,num_uniq])
        #num[i,3:num.shape[1]] = numpy.array( [ sum( [ (100.0 * float(g_n) / float(i+1)) >= c_f for g_n in genes ] ) for c_f in cfreq ] ) # core = found in at least x% of i
        #print(num)
        num = [num_all,num_new,num_uniq] + [ sum( [ (100.0 * float(g_n) / float(i+1)) >= c_f for g_n in genes ] ) for c_f in cfreq ]
        sys.stdout.write('%d\t%d\t%s\n' % (perm_i,i+1,'\t'.join( [str(x) for x in num] )))
    return perm_i, num

# NOTE MAIN
if __name__ == "__main__":
    # NOTE Args
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    # required
    parser.add_argument('--pam',    '-p', help='binary presence/absence matrix', required=True)
    parser.add_argument('--cfreq',  '-f', help='core gene min frequencies', nargs='+', default=[100,99,95,90])
    parser.add_argument('--ofile',  '-o', help='output file', default='/home/vgalata/tmp/tmp.txt')
    parser.add_argument('--perms',        help='number of permutations', default=10, type=int)
    parser.add_argument('--cores',        help='number of cores', default=1, type=int)
    parser.add_argument('--verbose','-v', help="Verbosity", action="store_true")
    args = parser.parse_args()
    
    # NOTE Init
    if args.verbose:
        sys.stdout.write('#Binary prs/abs matrix: %s\n#Core min. freq.: %s | %d permutations using %d cores\n' % ( args.pam,', '.join([str(x) for x in args.cfreq]),args.perms,args.cores ))
    
    # Read in the binary matrix
    pam = numpy.loadtxt(args.pam, delimiter='\t', skiprows=1, dtype='S')
    # remove first column (contains centroid IDs)
    pam = numpy.delete(pam, 0, 1)
    # convert to booleans
    pam = pam.astype(int)
    pam = pam.astype(bool)
    if args.verbose:
        sys.stdout.write('#Binary matrix: %d x %d\n' % (pam.shape[0],pam.shape[1]))
    #pam = pam[numpy.random.permutation(pam.shape[0])[0:10],1:10] # TEST
    #print(pam)
    
    sys.stdout.write('Perm\tN\tNumAll\tNumNew\tNumUniq\t%s\n' % '\t'.join([ 'NumCore%d' % c_f for c_f in args.cfreq ]))
    
    # Start permutations
    pool        = Pool(args.cores) # create pool for parallel computing
    pool_iter   = itertools.product(range(args.perms),[pam],[args.cfreq]) # iterable for the pool
    results     = pool.starmap( estimate_numbers , pool_iter ) # perform parallel computing
    pool.close(); pool.join() # wait until all finished and close the pool
    
    #with open(args.ofile,'w') as ofile:
        #ofile.write('Perm\tNumAll\tNumNew\tNumUniq\t%s\n' % '\t'.join([ 'NumCore%d' % c_f for c_f in args.cfreq ]))
        #for res in results:
           #for i in range(res[1].shape[0]):
               ##print( res[2][i,:] )
               #ofile.write( '%d\t%s\n' % ( res[0],'\t'.join( [ '%d' % x for x in res[1][i,:] ] ) ) )
    