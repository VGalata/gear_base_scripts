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
    parser$add_argument('-y_perf_files', help='Performance files', required=TRUE, nargs="+")
    parser$add_argument('-y_names', help='Names of phenotypes', required=TRUE, nargs="+")
    parser$add_argument('-y_models', help='model representation as text (e.g. feature list)', required=TRUE, nargs="+")
    parser$add_argument('-o_file', help='Output file', required=TRUE)
#     parser$add_argument('--src_path',   help='Print additional information.', default='.')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    
    return(parser)
}

## MAIN

args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
if(args$verbose){
    write(sprintf('Script summarizing CV results for all phenotypes\nGiven arguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}
# Check: same lengths
testit::assert(length(args$y_perf_files)==length(args$y_names))
testit::assert(length(args$y_perf_files)==length(args$y_models))

## Perf
perf_mat <- data.frame(Pheno=args$y_names, stringsAsFactors=FALSE, row.names=args$y_names, check.names=FALSE)

i <- 1
empty_perf_files <- c()
for(perf_file in args$y_perf_files){
    if(args$verbose) { write(sprintf("%s: %s",args$y_names[i],perf_file),stdout()) }
    # check name
    testit::assert(grepl(gsub('/| ','_',args$y_names[i]),perf_file))
    # if empty
    if(!file.exists(perf_file) | file.info(perf_file)$size==0){
        empty_perf_files <- c(empty_perf_files,i)
        i <- i+1
        next
    }
    p <- read.csv(file=perf_file, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    
    if(ncol(perf_mat)==1){
        perf_mat <- cbind(
            perf_mat,
            matrix(NA, nrow=nrow(perf_mat), ncol=ncol(p)-1, dimnames=list(rownames(perf_mat), colnames(p)[-1]))
        )
    }
    testit::assert( all( (colnames(p)[-1]) == (colnames(perf_mat)[-1]) ) )
    perf_mat[i,2:ncol(perf_mat)] <- p[p[,1]=='mean',-1]
    i <- i +1
}

# all files do not exists or are empty
if(length(empty_perf_files)==length(args$y_names)){
#     stop("All performance files were empty or do not exist")
    write("QUIT: All performance files were empty or do not exist", stdout())
    q("no", status=0, runLast=FALSE)
}

## Models
perf_mat$Model <- NA
i <- 1
empty_mod_files <- c()
for(mod_file in args$y_models){
    if(args$verbose) { write(sprintf("%s: %s",args$y_names[i],mod_file),stdout()) }
    # check name
    testit::assert(grepl(gsub('/| ','_',args$y_names[i]),mod_file))
    # if empty
    if(!file.exists(mod_file) | file.info(mod_file)$size==0){
        empty_mod_files <- c(empty_mod_files,i)
        i <- i+1
        next
    }
    mod <- sort(unique(readLines(mod_file, n=-1, skipNul=TRUE, warn=FALSE)))
    mod <- paste(mod, collapse=';')
    perf_mat[i,'Model'] <- mod
    i <- i+1
}

# all files do not exists or are empty
if(length(empty_mod_files)==length(args$y_names)){
#     stop("All model files were empty or do not exist")
    write("All model files were empty or do not exist", stdout())
    q("no", status=0, runLast=FALSE)
}

## Remove those with any empty file
empty_files <- unique(c(empty_perf_files,empty_mod_files))
if(length(empty_files)>0) {
    perf_mat <- perf_mat[-c(empty_files),,drop=FALSE]
    warning(sprintf('There are %d/%d/%d empty perf/mod/both files for %s',length(empty_perf_files), length(empty_mod_files), length(empty_files), paste(args$y_names[empty_files],collapse=' , ')))
}

## Remove columns with TN/TP/FN/FP
rm_cols <- c('TN','TP','FN','FP')
perf_mat <- perf_mat[,setdiff(colnames(perf_mat),rm_cols)]

## Sort by name
perf_mat <- perf_mat[sort(perf_mat$Pheno),]

## Save
write.table(x=perf_mat, file=args$o_file, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
