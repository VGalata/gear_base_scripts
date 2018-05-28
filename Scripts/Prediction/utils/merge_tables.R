#!/usr/bin/Rscript

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')
    # Required parameters: input/output files
    parser$add_argument('-mat_files', help='Matrix files', required=TRUE, nargs='+')
    parser$add_argument('-ofile', help='Output file', required=TRUE)
    # Samples
    parser$add_argument('--include_samples', help='')
    parser$add_argument('--exclude_samples', help='')
    # Other
    parser$add_argument('--src_path', help='', default='.')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    return(parser)
}

## MAIN
# Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
# Print info
if(args$verbose){
    write(sprintf('\nArguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

source(file.path(args$src_path,"my_utils.R"))

## Read in
mats <- list()
for(mat_file in args$mat_files){
    mat <- read.csv(file=mat_file, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
    mat <- t(mat)
    if(args$verbose){ write(sprintf("Read in: %s with %d variables and %d samples",mat_file,ncol(mat),nrow(mat)), stdout()) }
    mats<- c(mats, list(mat))
}

## Select samples
if(!is.null(args$include_samples) | !is.null(args$exclude_samples)){
    # read in sample lists
    include_samples <- exclude_samples <- NULL
    if(!is.null(args$include_samples)){ include_samples <- read_list(list_file=args$include_samples, verbose=args$verbose) }
    if(!is.null(args$exclude_samples)){ exclude_samples <- read_list(list_file=args$exclude_samples, verbose=args$verbose) }
    # select samples
    for(i in 1:length(mats)){
        mats[[i]] <- sel_mat_rows(mat=mats[[i]], incl_rnames=include_samples, excl_rnames=exclude_samples, verbose=args$verbose)
    }
    # check: all samples should be in each matrix
    if(!is.null(args$include_samples)){
        for(i in 1:length(mats)){
            testit::assert(sprintf(
                "Not all samples which should be included are there:\n%s\n",
                paste(include_samples[ sapply(include_samples, function(n){ !(n %in% rownames(mats[[i]])) }) ], collapse=',')
            ), all( sapply(include_samples, function(n){ n %in% rownames(mat) }) ))
        }
    }
}

## Samples should be in same order in each matrix
samples_intersect <- c()
for(i in 1:length(mats)){
    if(i==1){
        samples_intersect <- rownames(mats[[i]])
    }
    else{
        samples_intersect <- intersect(samples_intersect, rownames(mats[[i]]))
    }
}
if(args$verbose){ write(sprintf("Sample intersection: %d", length(samples_intersect)), stdout()) }
assert(length(samples_intersect)>0)

## Samples should be in same order in each matrix
samples_sorted <- sort(samples_intersect)
for(i in 1:length(mats)){
    testit::assert(all( sort(rownames(mats[[i]])) == samples_sorted ))
    mats[[i]] <- mats[[i]][samples_sorted,,drop=FALSE]
}

## Merge
big_mat <- NULL
for(mat in mats){
    if(is.null(big_mat)){
        big_mat <- mat
    }
    else{
        assert(length(intersect(colnames(big_mat),colnames(mat)))==0) # no intersection in column names
        assert(all(rownames(big_mat)==rownames(mat))) # same row names (assured by code above)
        big_mat <- cbind(big_mat, mat)
    }
}
if(args$verbose){
    write(sprintf("Final merged matrix has %d samples and %d variables.",nrow(big_mat),ncol(big_mat)), stdout())
}

## Save
write.table(x=t(big_mat), file=args$ofile, sep='\t', row.names=TRUE, col.names=TRUE, append=FALSE, quote=FALSE)