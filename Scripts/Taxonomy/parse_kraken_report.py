#!/usr/bin/python

# Modified version of cbackes script
__author__ = 'vgalata'

import argparse
import glob
import sys
import os
import time
from collections import OrderedDict
import itertools
from multiprocessing import Pool

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    # http://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python -> only change is "<" instead of "<="
    return abs(a-b) < max(rel_tol * max(abs(a), abs(b)), abs_tol)

def parse_kraken_report_file(filename, tax_rank='F', check_total=False, sampleID=None, orgmap=None):
    """
    :param filename: name of kraken report file
    FORMAT
    percentage_reads<tab>summed_up_reads<tab>direct_hits<tab>classification<tab>taxonomy_nr<tab>taxon_description<return>
    0.71    4236    4236    U   0   unclassified
    99.29   594353  3075    -   1   root
    98.64   590441  88  -   131567    cellular organisms
    98.59   590173  790 D   2       Bacteria
    98.07   587015  114 P   1224          Proteobacteria
    97.95   586336  157 C   1236            Gammaproteobacteria
    97.78   585320  0   O   91347             Enterobacteriales
    97.78   585320  73199   F   543             Enterobacteriaceae
    83.85   501902  8636    G   561               Escherichia
    82.16   491795  430394  S   562                 Escherichia coli
    2.47    14794   14794   -   362663                    Escherichia coli 536
      1) Percentage of reads covered by the clade rooted at this taxon
      2) Number of reads covered by the clade rooted at this taxon
      3) Number of reads assigned directly to this taxon
      4) A rank code, indicating (U)nclassified, (D)omain, (K)ingdom,
         (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
         All other ranks are simply '-'.
      5) NCBI taxonomy ID
      6) indented scientific name
    Reference: http://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports
    """
    if sampleID is None:
        try:
            sampleID = os.path.basename(filename).split("_")[0]
        except:
            sampleID = None
    
    rank_dict = {'D':'Domain','K':'Kingdom','P':'Phylum','C':'Class','O':'Order','F':'Family','G':'Genus','S':'Species'}
    
    best = dict.fromkeys(['taxon','pct','count'], None)
    expc = None
    if orgmap is not None and sampleID is not None and sampleID in orgmap:
        expc = {'taxon':orgmap[sampleID][0], 'rank':orgmap[sampleID][1], 'count':0, 'rank_count':0}
    rank_count = 0
    unclass_count = 0; unclass_pct = 0
    total_count = 0
    
    with open(filename) as infile:
        for line in infile:
            line = line.strip()
            if line:
                entries = line.split("\t")
                perc      = float(entries[0]) # Percentage of reads covered by the clade rooted at this taxon
                count     = int(entries[1]) # Number of reads covered by the clade rooted at this taxon
                count_dir = int(entries[2]) # Number of reads assigned directly to this taxon
                rank      = entries[3] # (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies; rest = '-'
                tax_name  = entries[5].strip() # NCBI taxon
                
                total_count += count_dir
                
                if rank == 'U': # unclassified
                    unclass_count = count
                    unclass_pct = perc
                    continue
                if expc is not None and rank in rank_dict and expc['rank'] == rank_dict[rank]: # expected taxon rank found
                    expc['rank_count'] += count
                    if tax_name == expc['taxon']: # expected taxon found
                        expc['count'] = count
                if rank == tax_rank: # needed rank found
                    rank_count += count # increase rank count
                    if best['taxon'] is None: # first
                        assert best['count'] is None and best['pct'] is None
                        best['taxon'] = tax_name
                        best['pct'] = perc
                        best['count'] = count
                    elif best['pct'] < perc: # found better
                        assert best['count'] < count, "Better pct but not count: %s, %f, %d vs. %s, %f, %d" % (best['taxon'], best['pct'], best['count'], tax_name, perc, count)
                        best['taxon'] = tax_name
                        best['pct'] = perc
                        best['count'] = count
                else:
                    pass
    
    if check_total:
        with open(filename) as infile:
            for line in infile:
                line = line.strip()
                if line:
                    entries = line.split("\t")
                    perc    = float(entries[0]) # Percentage of reads covered by the clade rooted at this taxon
                    count   = int(entries[1]) # Number of reads covered by the clade rooted at this taxon
                    perc2   = 100*float(count)/float(total_count) # compute own pct
                    assert isclose(perc, perc2, abs_tol=1e-1), "Total count = %d seems not to be correct: %.2f vs. %.2f" % (total_count, perc, perc2)
    
    if expc is None:
        expc = {'taxon':'NA', 'rank':'NA', 'count':0, 'rank_count':0}
    return best, rank_count, unclass_count, unclass_pct, total_count, expc, sampleID

def get_abundance_from_kraken_report(filename, tax_rank='F', min_abundance=0, check_total=False):
    rank_dict = {'D':'Domain','K':'Kingdom','P':'Phylum','C':'Class','O':'Order','F':'Family','G':'Genus','S':'Species'}
    
    taxa_dict = {}
    rank_count = 0
    unclass_count = 0; unclass_pct = 0
    total_count = 0
    
    with open(filename) as infile:
        for line in infile:
            line = line.strip()
            if line:
                entries = line.split("\t")
                perc      = float(entries[0]) # Percentage of reads covered by the clade rooted at this taxon
                count     = int(entries[1]) # Number of reads covered by the clade rooted at this taxon
                count_dir = int(entries[2]) # Number of reads assigned directly to this taxon
                rank      = entries[3] # (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies; rest = '-'
                tax_name  = entries[5].strip() # NCBI taxon
                
                total_count += count_dir
                
                if rank == 'U': # unclassified
                    unclass_count = count
                    unclass_pct = perc
                    continue
                if rank == tax_rank: # needed rank found
                    rank_count += count # increase rank count
                    if tax_name not in taxa_dict:
                        taxa_dict[tax_name] = [count, perc]
                    else:
                        print("WARNING: Taxon %s is already in dict and will be ignored")
                else:
                    pass
    # filter
    if min_abundance > 0:
        taxa = taxa_dict.keys() # all retrieved taxa
        for taxon in taxa:
            #taxon_pct = 100 * taxa_dict[taxon][0] / rank_count # percentage among all taxa of given rank
            #if taxon_pct < min_abundance: # below threshold
            if taxa_dict[taxon][1] < min_abundance: # below threshold
                del taxa_dict[taxon] # remove
    
    if check_total:
        with open(filename) as infile:
            for line in infile:
                line = line.strip()
                if line:
                    entries = line.split("\t")
                    perc    = float(entries[0]) # Percentage of reads covered by the clade rooted at this taxon
                    count   = int(entries[1]) # Number of reads covered by the clade rooted at this taxon
                    perc2   = 100*float(count)/float(total_count) # compute own pct
                    assert isclose(perc, perc2, abs_tol=1e-1), "Total count = %d seems not to be correct: %.2f vs. %.2f" % (total_count, perc, perc2)
    return taxa_dict, rank_count, unclass_count, unclass_pct, total_count

def read_orgmap(orgmap_file):
    orgmap = {} # empty dict
    with open(orgmap_file) as infile:
        for line in infile:
            line = line.strip()
            if line:
                entries = line.split("\t")
                if entries[0] not in orgmap:
                    orgmap[entries[0]] = entries[1:]
                else:
                    print("Sample %s already in dict!" % entries[0])
    return orgmap


if __name__ == "__main__":
    #parse command line arguments
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--report_folder', '-f', help='folder with kraken reports', required=True)
    parser.add_argument('--output', '-o', help='output file', required=True)
    parser.add_argument('--tax_rank', help='Fo which rank output should be created', default='S', choices=['D','K','P','C','O','F','G','S'])
    parser.add_argument("--check_total", help="whether computed total coun should be checked (makes things slower)", action="store_true")
    parser.add_argument('--orgmap_file', help='if given, pct and counts are retrieved only for given organisms correct organism for each sample: sample KielID, taxon, taxon level; tab. separated')
    parser.add_argument("--get_abundance", help="Instead of getting best taxon all taxa are retrieved with their counts", action="store_true")
    parser.add_argument("--min_abundance", help="If --get_abundance then taxa with at least this percentage among all taxa of given rank are retrieved", default=0, type=int)
    parser.add_argument('--cores', '-c', help='', default=1, type=int)
    args = parser.parse_args()
    
    report_folder = args.report_folder
    output_file = args.output
    
    if args.get_abundance and args.orgmap_file is not None:
        print("WARNING: Organism mapping file given and abundance should be trieved: Organism mapping file will be ignored")
    
    if not args.get_abundance:
        print("Extracting best assigned taxon at rank %s" % args.tax_rank)
        
        orgmap=None
        if args.orgmap_file is not None:
            orgmap = read_orgmap(args.orgmap_file)
        
        ordered_ranks = ['S', 'G', 'F', 'O', 'C', 'P', 'K', 'D','U']
        rank_dict = {'D':'Domain','K':'Kingdom','P':'Phylum','C':'Class','O':'Order','F':'Family','G':'Genus','S':'Species','U':'Unclassified'}
        rank_full = rank_dict[args.tax_rank]
        
        # OLD START
        #with open(output_file, 'w') as output:
            #sys.stdout.write("parsing: \n")
            #counter = 0
            #filelist = sorted(glob.glob(report_folder + "*.txt"))
            #for fname in filelist:
            ##for fname in filelist[0:3]: # TEST
                #sample_id = os.path.basename(fname).split("_")[0]
                ##print fname
                #sys.stdout.write("\r%d/%d: %s" % (counter, len(filelist),fname))
                #sys.stdout.write("\033[K") # Clear to the end of line
                #sys.stdout.flush()
                #time.sleep(1)
                #counter += 1
                #best, rank_count, unclass_count, unclass_pct, total_count, expc = parse_kraken_report_file(filename=fname, tax_rank=args.tax_rank, check_total=args.check_total, sampleID=sample_id, orgmap=orgmap)
                
                #if fname == filelist[0]: # first file -> add header
                    #output.write('ID\t%s\t%s_pct\t%s_count\t%s_rank_count\tUnclassified_pct\tUnclassified_count\tTotalCount\tExpectedTaxon\tExpectedTaxon_rank\tExpectedTaxon_count\tExpectedTaxon_rank_count\n' % (rank_full,rank_full,rank_full,rank_full))
                #output.write("%s\t%s\t%.2f\t%d\t%d\t%.2f\t%d\t%d\t%s\t%s\t%d\t%d\n" % (\
                    #sample_id, best['taxon'], best['pct'], best['count'], \
                    #rank_count, unclass_pct, unclass_count, total_count, \
                    #expc['taxon'], expc['rank'], expc['count'], expc['rank_count'])\
                #)
        # OLD END
        # NEW START
        reports     = sorted(glob.glob(report_folder + "*.txt"))
        pool        = Pool(args.cores) # create pool for parallel computing
        pool_iter   = itertools.product(reports,[args.tax_rank],[args.check_total],[None],[orgmap]) # filename, tax_rank='F', check_total=False, sampleID=None, orgmap=None
        results     = pool.starmap( parse_kraken_report_file , pool_iter ) # perform parallel computing
        pool.close(); pool.join() # wait until all finished and close the pool
        with open(output_file, 'w') as output:
            output.write('ID\t%s\t%s_pct\t%s_count\t%s_rank_count\tUnclassified_pct\tUnclassified_count\tTotalCount\tExpectedTaxon\tExpectedTaxon_rank\tExpectedTaxon_count\tExpectedTaxon_rank_count\n' % (rank_full,rank_full,rank_full,rank_full))
            for result in results:
                best, rank_count, unclass_count, unclass_pct, total_count, expc, sample_id  = result
                output.write("%s\t%s\t%.2f\t%d\t%d\t%.2f\t%d\t%d\t%s\t%s\t%d\t%d\n" % (\
                    sample_id, best['taxon'], best['pct'], best['count'], \
                    rank_count, unclass_pct, unclass_count, total_count, \
                    expc['taxon'], expc['rank'], expc['count'], expc['rank_count'])\
                )
        # NEW END
    else:
        print("Extracting all taxa at rank %s" % args.tax_rank)
        
        ordered_ranks = ['S', 'G', 'F', 'O', 'C', 'P', 'K', 'D','U']
        rank_dict = {'D':'Domain','K':'Kingdom','P':'Phylum','C':'Class','O':'Order','F':'Family','G':'Genus','S':'Species','U':'Unclassified'}
        rank_full = rank_dict[args.tax_rank]
        
        with open(output_file, 'w') as output:
            sys.stdout.write("parsing: \n")
            counter = 0
            filelist = sorted(glob.glob(report_folder + "*.txt"))
            for fname in filelist:
            #for fname in filelist[0:20]: # TEST
                sample_id = os.path.basename(fname).split("_")[0]
                taxa_dict, rank_count, unclass_count, unclass_pct, total_count = get_abundance_from_kraken_report(filename=fname, tax_rank=args.tax_rank, min_abundance=args.min_abundance, check_total=args.check_total)
                taxa_dict = OrderedDict(sorted(taxa_dict.items(), key=lambda t: -t[1][0])) # sort by count
                
                if fname == filelist[0]: # first file -> add header
                    output.write('ID\t%s\t%s_rank_count\tUnclassified_pct\tUnclassified_count\tTotalCount\n' % (rank_full,rank_full))
                output.write("%s\t%s\t%d\t%.2f\t%d\t%d\n" % (\
                    sample_id, \
                    ';'.join( ["%s:%d (%.2f)" % (taxon,taxon_c_p[0],taxon_c_p[1]) for taxon,taxon_c_p in taxa_dict.iteritems()] ), \
                    rank_count, unclass_pct, unclass_count, total_count \
                ))
                sys.stdout.write("\r%d/%d: %s" % (counter, len(filelist),fname))
                sys.stdout.write("\033[K") # Clear to the end of line
                sys.stdout.flush()
                time.sleep(1)
                counter += 1
    print("\nFinished\n")
