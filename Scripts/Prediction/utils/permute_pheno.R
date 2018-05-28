#!/usr/bin/Rscript

## Process binary feature file
## Assumption: values are 0/1, missing NA (but can be specified)

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('-pheno_file', help='', required=TRUE)
    parser$add_argument('-samples', help='', required=TRUE)
    parser$add_argument('-pheno_ofile', help='', required=TRUE)
    parser$add_argument('--src_path',   help='', default='.')
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

## LIBS
src_path <- args$src_path
source(sprintf("%s/utils/my_utils.R",src_path))

## READ IN
y <- read.csv(file=args$pheno_file, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)

samples <- read_list(list_file=args$samples, verbose=args$verbose)

testit::assert(all( sapply(samples, function(s){ s %in% rownames(y) }) ))

## Permute
y_perm <- y
for(i in 1:ncol(y_perm)){
    y_perm[samples,i] <- sample(x=y[samples,i])
}

## Save
write.table(x=y_perm, file=args$pheno_ofile, row.names=TRUE, col.names=TRUE, sep='\t', quote=FALSE)
