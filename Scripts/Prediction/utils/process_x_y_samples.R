#!/usr/bin/Rscript

## Process binary feature file
## Assumption: values are 0/1, missing NA (but can be specified)

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(gdsfmt))
suppressMessages(library(SNPRelate))

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('-geno_file', help='', required=TRUE)
    parser$add_argument('-pheno_file', help='', required=TRUE)
    parser$add_argument('-pheno_name', help='', required=TRUE)
    parser$add_argument('-o_samples', help='', required=TRUE)
    parser$add_argument('-o_stats', help='', required=TRUE)
    parser$add_argument('--include_samples', help='')
    parser$add_argument('--include_full', help='', action='store_true')
    parser$add_argument('--exclude_samples', help='')
    parser$add_argument('--rm_miss', help='', action='store_true')
    parser$add_argument('--miss_value', help='Value used for missing data', default="NA")
    parser$add_argument('--src_path', help='', default='.')
    parser$add_argument('--x_is_gds', help='', action='store_true')
    parser$add_argument('--verbose',  help='Print additional information.', action='store_true')
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
source(sprintf("%s/my_utils.R",src_path))

## READ IN
y <- read.csv(file=args$pheno_file, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
if(!args$x_is_gds){
    x <- read.csv(file=args$geno_file, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
    x <- t(x)
    if(args$miss_value!='NA'){
        x[x==args$miss_value]<-NA
        if(args$verbose){ write(sprintf("IMPORTANT: Missing value %s was replaced by NA.",args$miss_value),stdout()) }
    }
    # select samples
    if(!is.null(args$include_samples) | !is.null(args$exclude_samples)){
        include_samples <- exclude_samples <- NULL
        if(!is.null(args$include_samples)){ include_samples <- read_list(list_file=args$include_samples, verbose=args$verbose) }
        if(!is.null(args$exclude_samples)){ exclude_samples <- read_list(list_file=args$exclude_samples, verbose=args$verbose) }
        x <- sel_mat_rows(mat=x, incl_rnames=include_samples, excl_rnames=exclude_samples, verbose=args$verbose)
        y <- sel_mat_rows(mat=y, incl_rnames=include_samples, excl_rnames=exclude_samples, verbose=args$verbose)
        if(args$include_full){
            testit::assert(all( sapply(include_samples, function(s){ s %in% rownames(x) }) ))
            testit::assert(all( sapply(include_samples, function(s){ s %in% rownames(y) }) ))
        }
    }
}
if(args$x_is_gds){
    gds_file_open   <- snpgdsOpen(args$geno_file, readonly=TRUE)
    G_samples       <- read.gdsn(index.gdsn(gds_file_open, 'sample.id'))
    G_keep_samples  <- rep(TRUE, length(G_samples)) # samples to keep
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
                length(intersect(include_samples,exclude_samples)), sum(G_keep_samples)
            ), stdout())
        }
        if(sum(G_keep_samples)==0){ stop('There are no samples to select.') } # no samples remain
        y <- sel_mat_rows(mat=y, incl_rnames=include_samples, excl_rnames=exclude_samples, verbose=args$verbose)
        if(args$include_full){
            testit::assert(all( sapply(include_samples, function(s){ s %in% G_samples[G_keep_samples] }) ))
            testit::assert(all( sapply(include_samples, function(s){ s %in% rownames(y) }) ))
        }
    }
}

# remove samples with missing pheno
non_miss <- rep(TRUE, nrow(y))
if(args$rm_miss){
    testit::assert(args$pheno_name %in% colnames(y))
    non_miss <- !is.na(y[,args$pheno_name])
    write(sprintf("From %d samples %d have missing value for %s: %d remain",nrow(y),nrow(y)-sum(non_miss),args$pheno_name,sum(non_miss)), stdout())
}
# pheno stats
df <- data.frame(
    Class=c('miss','0','1'),
    Count=c(
        sum(is.na(y[non_miss,args$pheno_name])),
        sum(y[non_miss,args$pheno_name]==0,na.rm=TRUE),
        sum(y[non_miss,args$pheno_name]==1,na.rm=TRUE)
    )
)

if(args$verbose){
    write("Pheno stat.s:", stdout())
    print(df)
}

## OUTPUT
# Sample list
if(sum(non_miss)>0){
    write(x=sort(rownames(y)[non_miss]), file=args$o_samples, append=FALSE)
}
# Sample stats
write.table(x=df, file=args$o_stats, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
