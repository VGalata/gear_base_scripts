#!/usr/bin/Rscript

## Process binary feature file
## Assumption: values are 0/1, missing NA (but can be specified)

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(HiClimR))

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')
    # Required parameters: input/output files
    parser$add_argument('-mat_file', help='Binary feature file', required=TRUE)
    parser$add_argument('-o_features', help='', required=TRUE)
    parser$add_argument('-o_samples', help='')
    #
    parser$add_argument('--samples_header', help='samples in header', action='store_true')
    # Samples
    parser$add_argument('--include_samples', help='')
    parser$add_argument('--exclude_samples', help='')
    parser$add_argument('--include_features', help='')
    parser$add_argument('--exclude_features', help='')
    # Filtering
    parser$add_argument('--max_fq', help='Maximal (<=) frequency [%] of most frequent feature value, ignores missing.', type="double", default=1)
    parser$add_argument('--max_missfq', help='Maximal (>=) feature missing rate [%]', type="double", default=1)
    parser$add_argument('--miss_value', help='Value used for missing data', default="NA")
    # Correlation based filtering
    parser$add_argument('--rm_cor', help='Remove correlated features using correlation based hclust (d=1-cro^2)', action='store_true')
#     parser$add_argument('--rm_cor_method', help='', choices=c("pearson", "kendall", "spearman"), default='pearson')
    parser$add_argument('--rm_cor_hclust_method', help='', choices=c("complete","single","average","median","centroid","mcquitty"), default='average')
    parser$add_argument('--rm_cor_hclust_cut_h', help='At which height (d=1-cor^2) the hclust tree should be cut', type="double", default=1-(0.7**2))
    parser$add_argument('--rm_cor_hclust_cut_k', help='How many clusters should result after cutting hclust tree', type="integer", default=10)
    parser$add_argument('--rm_cor_hclust_use_h', help='Use height threshold to cut hclust tree', action='store_true')
    parser$add_argument('--rm_cor_hclust_cut_ofile', help='Output file to save corr. based clusters')
    # Other
    parser$add_argument('--src_path', help='', default='.')
    parser$add_argument('--cores', help='', default=2, type='integer')
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
mat <- read.csv(file=args$mat_file, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)

## Transpose if header contains samples
if(args$samples_header){
    mat <- t(mat)
}

## Init
# Check matrix values
testit::assert(all(sapply(sort(unique(as.vector(unlist(mat)))), function(x){ x %in% c(0,1,args$miss_value) })))

# set missing value to NA
if(args$miss_value!='NA'){
    mat[mat==args$miss_value]<-NA
    if(args$verbose){ write(sprintf("IMPORTANT: Missing value %s was replaced by NA.",args$miss_value),stdout()) }
}
print(mat[1:10,1:3])

## Select samples
if(!is.null(args$include_samples) | !is.null(args$exclude_samples)){
    include_samples <- exclude_samples <- NULL
    if(!is.null(args$include_samples)){ include_samples <- read_list(list_file=args$include_samples, verbose=args$verbose) }
    if(!is.null(args$exclude_samples)){ exclude_samples <- read_list(list_file=args$exclude_samples, verbose=args$verbose) }
    print(include_samples[1:10])
    mat <- sel_mat_rows(mat=mat, incl_rnames=include_samples, excl_rnames=exclude_samples, verbose=args$verbose)
    # all samples to include should be there
    if(!is.null(args$include_samples)){
        testit::assert(sprintf(
            "Not all samples which should be included are there:\n%s\n",
            paste(include_samples[ sapply(include_samples, function(n){ !(n %in% rownames(mat)) }) ], collapse=',')
        ), all( sapply(include_samples, function(n){ n %in% rownames(mat) }) ))
    }
    # if no sample output file given check if required
    if(!is.null(args$o_samples) & !all(rownames(mat)==setdiff(include_samples,exclude_samples))){
        stop(sprintf("File %s contains less/more samples than given sample list of samples to include",))
    }
}

## Select features
if(!is.null(args$include_features) | !is.null(args$exclude_features)){
    include_features <- exclude_features <- NULL
    if(!is.null(args$include_features)){ include_features <- read_list(list_file=args$include_features, verbose=args$verbose) }
    if(!is.null(args$exclude_features)){ exclude_features <- read_list(list_file=args$exclude_features, verbose=args$verbose) }
    mat <- sel_mat_cols(mat=mat, incl_cnames=include_features, excl_cnames=exclude_features, verbose=args$verbose)
}

## Filter columns
mat <- filter_mat_cols(mat=mat, missing_value=NA, max_fq=args$max_fq, max_missfq=args$max_missfq, cores=args$cores, verbose=args$verbose)

## Filter correlated
if(args$rm_cor){
    # Column (feature) correlations
#     mat_cor <- cor(mat, method=args$rm_cor_method)
    # http://finzi.psych.upenn.edu/library/HiClimR/html/fastCor.html --> pearson's corr. coeff.
    mat_cor <- fastCor(mat, nSplit=1, upperTri=FALSE, verbose=TRUE)
    
    # Plot
    if(!is.null(args$plot_pdffile)){ print(plot_col_cor(mat_cor)) }
    # Filter
    mat_cl <- remove_corr_hclust(
        mat=mat, mcor=mat_cor,
        cut_h=args$rm_cor_hclust_cut_h, cut_k=args$rm_cor_hclust_cut_k, use_h=args$rm_cor_hclust_use_h,
        aggl.method=args$rm_cor_hclust_method, verbose=args$verbose
    )
    mat <- mat_cl$mat
    mat_cl <- mat_cl$mat_cl
}

## Save
# Save table
# write.table(x=mat, file=args$o_file, sep='\t', quote=FALSE, append=FALSE, row.names=TRUE, col.names=TRUE)
# samples to keep
if(!is.null(args$o_samples)){
    write(x=rownames(mat),file=args$o_samples,  ncolumns=1, append=FALSE)
}
# features to keep
write(x=colnames(mat),file=args$o_features, ncolumns=1, append=FALSE)
# feature clusters if output file given
if (!is.null(args$rm_cor_hclust_cut_ofile)){
    write.table(x=mat_cl, file=args$rm_cor_hclust_cut_ofile, row.names=FALSE, sep='\t', quote=FALSE, append=FALSE)
}

if(args$verbose){
    write(sprintf('Finished\n'),stdout())
}
