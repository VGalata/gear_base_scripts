#!/usr/bin/python

__author__ = 'vgalata'

# NOTE IMPORTS
import sys
import argparse
import itertools
from multiprocessing import Pool
import logging
import pandas
from tqdm import tqdm

__HMM_TBLOUT_COLS__ = [
    'target_name', 'target_accession',
    'query_name', 'query_accession',
    'full_sequence_Evalue', 'full_sequence_score', 'full_sequence_bias',
    'best_1_domain_Evalue', 'best_1_domain_score', 'best_1_domain_bias',
    'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target'
]

__MAPS__ = None

def read_gene_tab(genes_tab):
    gene_map = pandas.read_csv(filepath_or_buffer=genes_tab, sep='\t', header=None, index_col=None)
    genes = {}
    centroids = {}
    for i in range(0, gene_map.shape[0]):
        c_id, c_genes = tuple(gene_map.iloc[i])
        c_genes = c_genes.split(';')
        genes.update(dict.fromkeys(c_genes, c_id))
        centroids[c_id] = c_genes
    logger.info('Centroid mapping for %d genes and %d centroids (%s)' % (len(genes), len(centroids), genes_tab))
    return {'centroids': centroids, 'genes': genes}

def read_hmm_hits():
    # read in
    hits = pandas.read_csv(filepath_or_buffer=args.hits, sep='\t', header=None, index_col=None)
    logger.info('Read in hits: %d entries' % hits.shape[0])
    # set column names
    hits.columns = __HMM_TBLOUT_COLS__
    # filter
    hits = hits.loc[ hits['rep'] > 0 ]
    logger.info('Filtered hits: %d remained' % hits.shape[0])
    return hits

def read_blast_hits():
    # read in
    hits = pandas.read_csv(filepath_or_buffer=args.hits, sep='\t', header=0, index_col=None)
    logger.info('Read in hits: %d entries' % hits.shape[0])
    # filter
    hits = hits.loc[ (hits['qcov'] >= args.min_cov) & (hits['scov'] >= args.min_cov) ]
    logger.info('Filtered hits: %d remained' % hits.shape[0])
    return hits

def get_centroid_gene_count(c_id):
    for mapping in __MAPS__:
        if c_id in mapping['centroids']:
            return len(mapping['centroids'][c_id])
    return None

# NOTE MAIN
if __name__ == "__main__":
    # Logger
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Args
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--genes', help='', required=True, nargs='+')
    parser.add_argument('--hits', help='', required=True)
    parser.add_argument('--ofile', help='', required=True)
    parser.add_argument('--min_cov', help='', type=int, default=0)
    parser.add_argument('--type', required=True, choices=['hmm', 'blast'])
    parser.add_argument('--s_id_col', choices=['query_name', 'query_accession'], default='query_name')
    parser.add_argument('--cores', help='', default=1, type=int)
    args = parser.parse_args()

    # Read in and filter hits
    hits = None
    if args.type == 'hmm':
        hits = read_hmm_hits()
    elif args.type == 'blast':
        hits = read_blast_hits()

    # Collect gene/centroid mappings
    pool        = Pool(args.cores)
    pool_iter   = itertools.product(args.genes)
    __MAPS__    = pool.starmap( read_gene_tab , pool_iter )
    pool.close(); pool.join()
    logger.info('Collected gene-centroid mappings.')

    # Count
    centroid_subject_hits = {} # Centroid, subject -> #hits
    for i in tqdm(range(0, hits.shape[0])):
        # gene and subject id
        g_id = None; s_id = None
        if args.type == 'hmm':
            g_id = hits.iloc[i]['target_name']
            s_id = hits.iloc[i][args.s_id_col]
        elif args.type == 'blast':
            g_id = hits.iloc[i]['qseqid']
            s_id = hits.iloc[i]['sseqid']
        # find mapping
        for mapping in __MAPS__:
            # mapping found
            if g_id in mapping['genes']:
                # centroid id
                c_id = mapping['genes'][g_id]
                # add to dict
                if (c_id, s_id) not in centroid_subject_hits:
                    centroid_subject_hits[(c_id,s_id)] = 0
                centroid_subject_hits[(c_id,s_id)] = centroid_subject_hits[(c_id,s_id)] + 1
                break

    with open(args.ofile, 'w') as ofile:
        ofile.write('CentroidID\tSubjectID\tCount\tCentroidGeneCount\n')
        for key, count in centroid_subject_hits.items():
            ofile.write('%s\t%s\t%d\t%d\n' % (key[0], key[1], count, get_centroid_gene_count(key[0])))
