#!/usr/bin/Rscript

## Process binary feature file
## Assumption: values are 0/1, missing NA (but can be specified)

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
source("my_utils.R")
# suppressMessages(library(ggplot2))
# source("/home/vgalata/Bacteria_pipeline/Analysis_utils/ggplot2_multiplot.R")

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
	parser <- ArgumentParser()
	# Required parameters: input/output files
	parser$add_argument('-mat_file', help='Binary phenotype file', required=TRUE)
	parser$add_argument('-o_file', help='Output file', required=TRUE)
	# Samples
	parser$add_argument('--include_samples', help='Sample list (one sample per line, comment char is hash/diamond)')
	parser$add_argument('--exclude_samples', help='Sample list (one sample per line, comment char is hash/diamond)')
	# Filtering
	parser$add_argument('--max_fq', help='Maximal (<=) frequency [%] of most frequent feature value.', type="double", default=1)
	parser$add_argument('--max_missfq', help='Remove SNPs missing in at least (>=) x percent [%] of samples.', type="double", default=1)
	parser$add_argument('--missing_value', help='Value used for missing data', default="NA")
	# Other
	parser$add_argument('--plot_pdffile', help='If given plots are generated as saved in the file')
	parser$add_argument('--verbose', help='Print additional information.', action='store_true')
	
	return(parser)
}

## MAIN
## Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
if(args$missing_value=='NA'){ args$missing_value <- NA }
# Print info
if(args$verbose){
	write(sprintf('Script for processing binary phenotype data\nGiven arguments:'),stdout())
	write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

## Init
mat <- read_y(y_file=args$mat_file, verbose=args$verbose)
# Check values
mat_u_values <- sort(unique(as.vector(unlist(mat)))) # unique values in matrix
allowed_values <- c(0,1,args$missing_value) # allowed values
testit::assert(all(sapply(mat_u_values, function(x){ x %in% allowed_values })))
# set missing value to NA
if(!is.na(args$missing_value)){
	mat[mat==args$missing_value]<-NA
	if(args$verbose){ write(sprintf("IMPORTANT: Missing value %s was replaced by NA.",args$missing_value),stdout()) }
}
# Plot pdf
if(!is.null(args$plot_pdffile)){ pdf(file=args$plot_pdffile) }

## Select samples
if(!is.null(args$include_samples) | !is.null(args$exclude_samples)){
	# Samples to iclude/exclude
	include_samples <- exclude_samples <- NULL
	if(!is.null(args$include_samples)){ include_samples <- read_list(list_file=args$include_samples, verbose=args$verbose) }
	if(!is.null(args$exclude_samples)){ exclude_samples <- read_list(list_file=args$exclude_samples, verbose=args$verbose) }
	# Select
	mat <- sel_mat_rows(mat=mat, incl_rnames=include_samples, excl_rnames=exclude_samples, verbose=args$verbose)
}

## Filter samples/rows
## Not implemeneted because not needed now (no missing data for some samples, only for all)

## Filter columns
# Plot
if(!is.null(args$plot_pdffile)){ print(plot_phenos(mat)) }
# Filter
mat <- filter_mat_cols(mat=mat, missing_value=NA, max_fq=args$max_fq, max_missfq=args$max_missfq, verbose=args$verbose)
# Plot
if(!is.null(args$plot_pdffile)){ print(plot_phenos(mat)) }

## Save
# Clode PDF
if(!is.null(args$plot_pdffile)){ dev.off() }
# Write table
write.table(x=mat, file=args$o_file, sep='\t', quote=FALSE, append=FALSE, row.names=TRUE, col.names=TRUE)