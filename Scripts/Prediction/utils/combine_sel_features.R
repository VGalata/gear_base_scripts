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
    parser$add_argument('-sel_features_files', help='', required=TRUE, nargs="+")
    parser$add_argument('-o_file', help='Output file', required=TRUE)
    # other
    parser$add_argument('--min_rate', help='', type='double', default=90)
    parser$add_argument('--src_path',   help='Print additional information.', default='.')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    
    return(parser)
}


## MAIN

## Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
if(args$verbose){
    write(sprintf('Combine selected features\nGiven arguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

source(sprintf("%s/my_utils.R",args$src_path))

## Do work
# Read in lists
sel_features_list <- list()
for(sel_features_file in args$sel_features_files){
    if(!file.exists(sel_features_file) | file.info(sel_features_file)$size==0){
        write(sprintf('Could not read features from %s', args$features))
        next
    }
    sel_features <- read_list(list_file=sel_features_file, verbose=args$verbose)
    if(is.null(sel_features)){
        write(sprintf('MISS: Could not read features from %s', args$features))
        next
    }
    if(args$verbose){ write(sprintf("File %s: %d features", sel_features_file, length(sel_features)), stdout()) }
    
    sel_features_list <- c(sel_features_list, list(sel_features))
    names(sel_features_list)[length(sel_features_list)] <- sel_features_file
}

if(length(sel_features_list)==0){
    write("QUIT: All feature lists do not exist or are empty", stdout())
    q("no", status=0, runLast=FALSE)
}

# All features appearing in any list
sel_features <- sort(unique(unlist(sel_features_list)))
# Create table
sel_features_tab <- data.frame(
    matrix(NA, nrow=length(sel_features), ncol=length(sel_features_list), dimnames=list(sel_features,names(sel_features_list))),
    row.names=sel_features,
    stringsAsFactors=FALSE, check.names=FALSE
)

for(f in colnames(sel_features_tab)){
    sel_features_tab[sel_features_list[[f]],f] <- 1
}

# Select
# without NA entry <=> contained in file
o_features <- rownames(sel_features_tab)[apply(sel_features_tab,1,function(x){ (100*sum(!is.na(x))/length(x)) >= args$min_rate })]
if(args$verbose){ write(sprintf("From %d features %d were selected", nrow(sel_features_tab), length(o_features)), stdout()) }

## Save
write(x=o_features, file=args$o_file)
write.table(sel_features_tab, file=sprintf("%s.tab",args$o_file), row.names=TRUE, col.names=TRUE, sep='\t', quote=FALSE)