#!/usr/bin/Rscript

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(tools))



## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    # required parameters
    parser$add_argument('-perf_files', help='CV performance files', required=TRUE, nargs="+")
    parser$add_argument('-o_file', help='Output file', required=TRUE)
    parser$add_argument('--src_path',   help='Print additional information.', default='.')
    parser$add_argument('--reps',   help='Print additional information.', required=TRUE, type='integer')
    parser$add_argument('--min_rep_pct',   help='Print additional information.', default=50.0, type='double')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    
    return(parser)
}

## MAIN

## Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
if(args$verbose){
    write(sprintf('Script summarizing CV prediction results\nGiven arguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

source(sprintf("%s/perf_utils.R",args$src_path))

## Read
perf_mat <- NULL
for(perf_file in args$perf_files){
	if(!file.exists(perf_file) | file.info(perf_file)$size==0){
# 		stop(sprintf("Empty performance file %s",perf_file))
        write(sprintf("MISS: Empty performance file %s",perf_file))
        next
	}
    p <- read.csv(file=perf_file, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    if(args$verbose){
        write(sprintf("Performance file %s: %d rows, %d cols",perf_file,nrow(p),ncol(p)),stdout())
    }
    if(is.null(perf_mat)){
        perf_mat <- p
    }
    else{
        testit::assert(all(colnames(p)==colnames(perf_mat)))
        perf_mat <- rbind(perf_mat,p)
    }
    rownames(perf_mat)[nrow(perf_mat)] <- regmatches(perf_file, regexpr('rep[[:digit:]]+',perf_file))
}

## Check whether there are too few repetitions
# number of reps in table
reps <- sum(grepl('^rep[0-9]+',rownames(perf_mat)))
if((100*reps/args$reps) < args$min_rep_pct){
#     stop(sprintf("No performance because only %d repetitions have results",reps))
    write(sprintf("QUIT: No performance because only %d repetitions have results",reps), stdout())
    q("no", status=0, runLast=FALSE)
}

## Check: same number of samples in both classes
testit::assert(all(perf_mat$num0==mean(perf_mat$num0)))
testit::assert(all(perf_mat$num1==mean(perf_mat$num1)))

## Mean, sd
perf_mat <- rbind(
    perf_mat,
    mean=apply(perf_mat, 2, mean),
    sd=apply(perf_mat, 2, sd)
)

## First column
perf_mat <- data.frame(
    ID=rownames(perf_mat), perf_mat,
    check.names=FALSE, stringsAsFactors=FALSE
)

## Save
write.table(x=perf_mat, file=args$o_file, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)