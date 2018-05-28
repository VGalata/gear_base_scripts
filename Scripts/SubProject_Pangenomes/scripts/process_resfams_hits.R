#!/usr/bin/Rscript

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))

## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--hits', help='', required=TRUE)
    parser$add_argument('--annot', help='', required=TRUE)
    parser$add_argument('--meta_rds', help='', required=TRUE, nargs='+')
    parser$add_argument('--obname', help='', required=TRUE)
    parser$add_argument('--src_path', help='')
    return(parser)
}

## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## Help functions
suppressMessages(source(sprintf('%s/utils.R',args$src_path)))

## Log file
log_file <- sprintf('%s.log', args$obname)
init_log(log_file, args)

## Data
# sample gene - essential gene hits
hits <- read.csv(file=args$hits, header=FALSE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
colnames(hits) <- hmm_tblout_cols # set column names
hits$sample <- sapply(hits[,'target_name'], geneID2sampleID) # add sample ID (derived from target ID)
# resfams annotations
annot <- read.csv(file=args$annot, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
rownames(annot) <- annot[,'Resfam ID']
# sample meta-data
tmp <- read_meta_rds(args$meta_rds)
meta <- tmp$meta
samples <- tmp$samples
write(sprintf('Meta data: %d samples (%d, rest %d)', nrow(meta), length(samples$main), length(samples$rest)), log_file, append=TRUE)

## Main taxa
# From sampels used in pan-genomes get unique species taxa
main_taxa <- convert_taxa(sort(unique(meta[samples$main, 'TaxSpecies'])))
write(sprintf('Main species taxa: %s', paste(main_taxa, collapse=', ')), log_file, append=TRUE)

## Summarize
# Pre-filter hits: at least one reported domain
# hits <- hits[hits$rep > 0,]
# write(sprintf('Filtered by \"rep\" > 0: %d hits remained', nrow(hits)), log_file, append=TRUE)
# Reformat hits to data.table format
hits_dt <- data.table(hits)
# Number of hits per query and sample
hits_sum <- aggregate(rep(1,nrow(hits)), by=list(query_name=hits[,'query_accession'], sample=hits$sample), FUN=sum)
write('Aggregated hit counts', log_file, append=TRUE)

# Table: Samples x genes = number of hits
# reshape: sample x query
ov_g_s <- dcast(hits_sum, sample ~ query_name)
# set sample IDs as row names and remove as column
rownames(ov_g_s) <- ov_g_s$sample
ov_g_s <- ov_g_s[,colnames(ov_g_s)!='sample']
# add missing resfam genes
for(q_n in rownames(annot)){
    if(! q_n %in% colnames(ov_g_s)){
        ov_g_s[,q_n] <- 0
    }
}
# add missing samples (e.g. no/bad assembly)
for(s_id in rownames(meta)){
    if(! s_id %in% rownames(ov_g_s)){
        ov_g_s[s_id,] <- 0
    }
}
# replace NAs by 0 (NAs produced by dcast = no hit)
ov_g_s[is.na(ov_g_s)] <- 0
# checks
testit::assert(nrow(ov_g_s) == nrow(meta))
testit::assert(ncol(ov_g_s) == nrow(annot))
# NOTE keep only samples used in pan-genomes !!!
ov_g_s <- ov_g_s[samples$main,]
write(sprintf('\nKeeping only relevant samples: %d', nrow(ov_g_s)), log_file, append=TRUE)
# presence pct. for each e. gene
ov_g_s_prs <- apply(ov_g_s, 2, function(x){ 100 * sum(x > 0) / nrow(ov_g_s) })

# Table: sample, taxon, number of total/unique/mult. hits; total = unique + mult.
ov_n_s <- data.frame(
    sample=rownames(ov_g_s),
    taxon=convert_taxa(meta[rownames(ov_g_s), 'TaxSpecies']),
    total=apply(ov_g_s, 1, function(x){ sum(x > 0) }),
    uniq= apply(ov_g_s, 1, function(x){ sum(x == 1) }),
    mult= apply(ov_g_s, 1, function(x){ sum(x >= 2) }),
    row.names=rownames(ov_g_s), stringsAsFactors=FALSE, check.names=FALSE
)
# discard not relevant taxa
ov_n_s$taxon <- sapply(as.character(ov_n_s$taxon), function(x){ ifelse(x %in% main_taxa, x, 'Other') })
# reshape: sample ID, taxon, stat.s type (total/unique/mult), number of hits for stat.s type
ov_n_s <- reshape::melt.data.frame(ov_n_s, id.vars=c('sample', 'taxon'))

# Table: species x resfams = mean
# add sample ID and taxon, and reshape -> sample ID, taxon, query name, XX
ov_g_sp <- reshape::melt.data.frame(
    cbind(
        ov_g_s,
        sample=rownames(ov_g_s), # sample ID
        taxon=convert_taxa(meta[rownames(ov_g_s), 'TaxSpecies']) # taxon
    ),
    id.vars=c('sample', 'taxon')
)
# discard not relevant taxa
ov_g_sp$taxon <- sapply(as.character(ov_g_sp$taxon), function(x){ ifelse(x %in% main_taxa, x, 'Other') })
# aggregate by taxon and XX: compute mean of XX
ov_g_sp <- aggregate(ov_g_sp$value, by=list(taxon=ov_g_sp$taxon, query_name=ov_g_sp$variable), FUN=mean)
# reshape: taxon x query
ov_g_sp <- dcast(ov_g_sp, taxon ~ query_name)
# set row names and remove sample IDs as column
rownames(ov_g_sp) <- ov_g_sp$taxon
ov_g_sp <- ov_g_sp[,colnames(ov_g_sp) != 'taxon']

# Save
objects <- list(
    ov_g_s=ov_g_s,
    ov_g_sp=ov_g_sp,
    ov_n_s=ov_n_s,
    ov_g_s_prs=ov_g_s_prs
)
saveRDS(objects, sprintf('%s.rds', args$obname))
