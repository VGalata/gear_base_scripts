#!/usr/bin/Rscript

## Process GDS file: filter by samples and SNPs
## Assumption: genotypes are 0/1 (ALT/REF), missing > 2

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(gdsfmt))
suppressMessages(library(SNPRelate))
suppressMessages(library(ggplot2))

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    # Required parameters
    parser$add_argument('-gds_file', help='GDS file (path/name)', required=TRUE)
    parser$add_argument('-o_features', help='', required=TRUE)
    parser$add_argument('-o_samples', help='')
    # Samples
    parser$add_argument('--include_samples', help='Sample list (one sample per line, comment char is hash/diamond)')
    parser$add_argument('--exclude_samples', help='Sample list (one sample per line, comment char is hash/diamond)')
    # SNP Filtering
    parser$add_argument('--max_fq', help='Maximal frequency of most frequent feature value.', type="double", default=1)
    parser$add_argument('--max_missfq', help='Maximal feature missing rate/frequency.', type="double", default=1)
    # LD pruning
    parser$add_argument('--ld_prun', help='Whether LD pruning should be done', action='store_true')
    parser$add_argument('--ld_prun_method', help='', choices=c("composite", "r", "dprime", "corr"), default='corr')
    parser$add_argument('--ld_prun_slide_max_bp', help='', type='integer', default=500000)
    parser$add_argument('--ld_prun_slide_max_n', help='Set to -1 (default) if it should be NA.', type='integer', default=-1)
    parser$add_argument('--ld_prun_threshold', help='', type='double', default=0.2)
    # other
    parser$add_argument('--src_path', help='', default='.')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    
    return(parser)
}

## MAIN
## Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
if(args$ld_prun_slide_max_n==-1){ args$ld_prun_slide_max_n <- NA }
# Print info
if(args$verbose){
    write(sprintf('Script for processing GDS objects containing VCF data\nGiven arguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

source(file.path(args$src_path,"my_utils.R"))

## Init
# Open file
gds_file_open <- snpgdsOpen(args$gds_file, readonly=TRUE)
# Genotype matrix, samples and snp IDs
G         <- index.gdsn(gds_file_open, 'genotype')
G_samples <- read.gdsn(index.gdsn(gds_file_open, 'sample.id'))
G_snps    <- read.gdsn(index.gdsn(gds_file_open, 'snp.id'))
# Check genotype values
G_u_values <- sort(unique(as.vector(read.gdsn(G)))) # unique values in the matrix (sorted)
allowed_values <- c(0,1,2,3)
testit::assert(
    sprintf("GDS genotype matrix: Values %s are not all in allowed values %s.",paste(G_u_values,collapse=', '),paste(allowed_values,collapse=', ')),
    all(sapply(G_u_values, function(x){ x %in%  allowed_values }))
)

## Processing
## Select samples from given list
G_keep_samples <- rep(TRUE, length(G_samples)) # samples to keep
if(!is.null(args$include_samples) | !is.null(args$exclude_samples)){
    # Sample lists
    if(!is.null(args$include_samples)){ include_samples <- read_list(list_file=args$include_samples, verbose=args$verbose) }
    else{ include_samples <- G_samples } # include all
    if(!is.null(args$exclude_samples)){ exclude_samples <- read_list(list_file=args$exclude_samples, verbose=args$verbose) }
    else{ exclude_samples <- c() } # exclude none
    # Select
    G_keep_samples  <- sapply(G_samples, function(x){ (x %in% include_samples) & (!x %in% exclude_samples) })
    if(args$verbose){
        write(sprintf(
            '\nSelecting samples: GDS contains %d samples\n\t%d to include, %d are in GDS\n\t%d to exclude, %d are in GDS\n\t%d include and exclude\n\t%d remained in GDS',
                length(G_samples),
                length(include_samples), length(intersect(G_samples,include_samples)),
                length(exclude_samples), length(intersect(G_samples,exclude_samples)),
                length(intersect(include_samples,exclude_samples)),
                sum(G_keep_samples)
        ), stdout())
    }
    if(sum(G_keep_samples)==0){ stop('There are no samples to select.') } # no samples remain
    # all to include are there
    if(!is.null(args$include_samples)){
        testit::assert(sprintf(
            "Not all given samples are in matrix after processing in/ex samples:\n%s\n",
            paste(include_samples[ sapply(include_samples, function(n){ !(n %in% G_samples[G_keep_samples]) }) ], collapse=',')
        ),all( sapply(include_samples, function(n){ n %in% G_samples[G_keep_samples] }) ))
    }
}

## Filter SNPs: use selected samples and _all_ SNPs (selection arg.)
G_keep_snps <- apply.gdsn(
    node=G, margin=2, FUN=function(x){
        if(all(x>2)){ return(FALSE) } # all missing -> do not keep
        x_mm <- length(unique(x[x<=2]))==1 # TRUE if monomorphic (only one non-missing value)
        x_tab<- table(x[x<=2]) # table without missing
        x_mf <- 100*max(x_tab)/sum(x_tab) # maximal genotype freq. [%]
        x_mr <- 100*sum(x>2)/length(x) # missing rate [%]
        # keep if: non-monomorphic AND max. freq. <= t1 AND missing rate <= t2
        return((!x_mm) & x_mf<=args$max_fq & x_mr<=args$max_missfq)
    },
    selection=list(G_keep_samples,rep(TRUE, length(G_snps))), as.is='logical'
)
if(args$verbose){
    write(sprintf(
        "\nGDS: SNP filtering: From %d features, %d were removed (monomorphic, max. freq. rate > %0.3f, missing rate > %0.3f): %d remained",
        length(G_snps), length(G_snps)-sum(G_keep_snps), args$max_fq, args$max_missfq, sum(G_keep_snps)
    ),stdout())
}

## LD pruning
if(args$ld_prun){
    # LD pruning
    G_keep_ld_snps <- unlist(snpgdsLDpruning(
        gdsobj=gds_file_open, sample.id=G_samples[G_keep_samples], snp.id=G_snps[G_keep_snps], autosome.only=FALSE,
        remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, # should not remove anything as already done above
        method=args$ld_prun_method, slide.max.bp=args$ld_prun_slide_max_bp, slide.max.n=args$ld_prun_slide_max_n, ld.threshold=args$ld_prun_threshold,
        num.thread=1, # parallel version is not implemented
        verbose=FALSE #args$verbose set to FALSE because otherwise prints for each chrom. (i.e. many output for pan-genomes)
    ))
    # Check whether all pruned SNPs are in selected SNPs
    testit::assert("(pruned SNPs) != Intersect(selected SNPs, pruned SNPs)", length(G_keep_ld_snps)==length(intersect(G_snps[G_keep_snps],G_keep_ld_snps)))
    # For 1st to p-th SNP:  If keep-flag then (check whether ID in pruned) else FALSE
    G_keep_snps <- sapply(1:length(G_snps), function(x){ ifelse(!G_keep_snps[x],FALSE,G_snps[x] %in% G_keep_ld_snps) })
    # Checks
    testit::assert( length(G_keep_ld_snps) == sum(G_keep_snps) ) # same number
    testit::assert( sort(G_keep_ld_snps)   == sort(G_snps[G_keep_snps]) ) # same SNPs
    if(args$verbose){
        write(sprintf( "\nGDS: LD prunning: %d features remained", sum(G_keep_snps) ),stdout())
    }
}

## Subselection of G w.r.t. selected samples and SNPs
# # Subset of matrix by samples and SNPs
# G_sub <- readex.gdsn(G, list(G_keep_samples, G_keep_snps))
# # Set dim. names
# dimnames(G_sub) <- list(G_samples[G_keep_samples], G_snps[G_keep_snps])
# # Missing -> NA
# G_sub[G_sub>2] <- NA # set missing values to NAs
# if(args$verbose){ write('IMPORTANT: Missing value (genotype > 2, usually 3) was replaced ba NA',stdout()) }

## Load gene presence/absence

## Save and close
# # Save matrix
# write.table(x=G_sub, file=args$o_file, sep='\t', quote=FALSE, append=FALSE, row.names=TRUE, col.names=TRUE)
# samples to keep
if(!is.null(args$o_samples)){
    write(x=G_samples[G_keep_samples], file=args$o_samples,  ncolumns=1, append=FALSE)
}
# features to keep
write(x=G_snps[G_keep_snps], file=args$o_features, ncolumns=1, append=FALSE)
# Close GDS file
snpgdsClose(gds_file_open)