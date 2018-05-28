#!/usr/bin/Rscript

## Convert a VCF file to a genotype matrix stored in GDS format:
##  biallelic.only: ?
##  copy.num.of.ref: 0 = Alternative, 1 = Reference (haploid); 0,1,2 for diploid
## Can rename samples if Kiel seq. sample ID was used

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(gdsfmt))
suppressMessages(library(SNPRelate))

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('-vcf_file', help='VCF file which should be converted (path/name)', required=TRUE)
    parser$add_argument('-gds_file', help='Output GDS file (path/name)', required=TRUE)
    parser$add_argument('-bin_file', help='Output GDS file (path/name) for bin. representation', required=TRUE)
    parser$add_argument('--include_samples', help='Sample list (one sample per line, comment char is hash/diamond)')
    parser$add_argument('--exclude_samples', help='Sample list (one sample per line, comment char is hash/diamond)')
    parser$add_argument('--snpgdsVCF2GDS_method', choices=c("biallelic.only", "copy.num.of.ref"), default="copy.num.of.ref")
    parser$add_argument('--no_rename', help='Whether sample names should be kept as they are', action='store_true')
    parser$add_argument('--src_path', help='', default='.')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    return(parser)
}

## MAIN
## Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
# Print info
if(args$verbose){
    write(sprintf('Script for converting VCF to GDS (Bioconductor \"gdsfmt\" for storing genotype data)\nGiven arguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

source(file.path(args$src_path,"my_utils.R"))

## Convert
if(args$verbose){ write(sprintf("Starting VCF -> GDS: %s",format(Sys.time(), "%Y.%m.%d %X")), stdout()) }

snpgdsVCF2GDS(vcf.fn=args$vcf_file, out.fn=args$gds_file, method=args$snpgdsVCF2GDS_method, verbose=args$verbose, snpfirstdim=FALSE)

if(args$verbose){ write(sprintf("Finished VCF -> GDS: %s",format(Sys.time(), "%Y.%m.%d %X")), stdout()) }

## Open for further processing
gds_file_open   <- snpgdsOpen(args$gds_file, readonly=FALSE)

## Replace sample names
if(!args$no_rename){
    gds_sampleIDs   <- index.gdsn(gds_file_open, index='sample.id') # get sample IDs from GDS
    sampleIDs       <- read.gdsn(gds_sampleIDs) # sample IDs as vector
    # New sample IDs: Example: "C3832-407_GGACTCCT-TAGATCGC_L004_R" -> "C3832-407"
#   new_sampleIDs <- gsub(pattern='_(A|T|C|G|U)+-(A|T|C|G|U)+_L[[:digit:]]+_R.*', replacement='', x=sampleIDs) ## may need to adjust that
    new_sampleIDs   <- sub(pattern='_R$', replacement='', x=sampleIDs) ## may need to adjust that
    rename.gdsn(gds_sampleIDs, "old.sample.id") # Rename atrribute in GDS containing sample IDs
    gds_new_sampleIDs <- add.gdsn(gds_file_open, "sample.id", val=new_sampleIDs, storage='string') # Add new sample IDs, attribute name is the name of the former sample ID attribute (wich was renamed)
    # Print info
    if(args$verbose){
        write(sprintf(
            "Samples were renamed:\n\tFrom %s ...\n\tTo %s ...\n",
            paste(read.gdsn(index.gdsn(gds_file_open,'old.sample.id'), start=1, count=5),collapse=', '),
            paste(read.gdsn(index.gdsn(gds_file_open,'sample.id'), start=1, count=5),collapse=', ')
        ), stdout())
    }
}

## Set variant (SNP) names: Chrom. + pos. made unique
## Not used because ID is then too long for EIGENSTRAT
# gds_snpIDs      <- index.gdsn(gds_file_open, index='snp.id') # get sample IDs from GDS
# snpIDs          <- read.gdsn(gds_snpIDs) # sample IDs as vector
# snpChrom        <- read.gdsn( index.gdsn(gds_file_open, index='snp.chromosome') )
# snpPos          <- read.gdsn( index.gdsn(gds_file_open, index='snp.position') )
# snpAllele       <- read.gdsn( index.gdsn(gds_file_open, index='snp.allele') )
# new_snpIDs      <- make.unique(paste(snpChrom,snpPos,gsub('/|,','_',snpAllele), sep='_'), sep='_')
# rename.gdsn(gds_snpIDs, "old.snp.id") # Rename atrribute in GDS containing sample IDs
# gds_new_snpIDs  <- add.gdsn(gds_file_open, "snp.id", val=new_snpIDs, storage='string') # Add new sample IDs, attribute name is the name of the former sample ID attribute (wich was renamed)
# # Print info
# if(args$verbose){
# write(sprintf(
#     "Variants were renamed:\n\tFrom %s ...\n\tTo %s ...\n",
#         paste(read.gdsn(index.gdsn(gds_file_open,'old.snp.id'), start=1, count=5),collapse=', '),
#         paste(read.gdsn(index.gdsn(gds_file_open,'snp.id'), start=1, count=5),collapse=', ')
#     ), stdout())
# }

## Set variant (SNP) names: rs[i], i = 1 .. n
gds_snpIDs      <- index.gdsn(gds_file_open, index='snp.id') # get sample IDs from GDS
snpIDs          <- read.gdsn(gds_snpIDs) # sample IDs as vector
new_snpIDs      <- paste('var_',snpIDs,sep='')
rename.gdsn(gds_snpIDs, "old.snp.id") # Rename atrribute in GDS containing sample IDs
gds_new_snpIDs  <- add.gdsn(gds_file_open, "snp.id", val=new_snpIDs, storage='string') # Add new sample IDs, attribute name is the name of the former sample ID attribute (wich was renamed)
# Print info
if(args$verbose){
write(sprintf(
    "Variants were renamed:\n\tFrom %s ...\n\tTo %s ...\n",
        paste(read.gdsn(index.gdsn(gds_file_open,'old.snp.id'), start=1, count=5),collapse=', '),
        paste(read.gdsn(index.gdsn(gds_file_open,'snp.id'), start=1, count=5),collapse=', ')
    ), stdout())
}

## Save GDS as binary matrix w.r.t. given samples
# Genotype matrix, samples, SNPs
G         <- index.gdsn(gds_file_open, 'genotype')
G_samples <- read.gdsn(index.gdsn(gds_file_open, 'sample.id'))
G_snps    <- read.gdsn(index.gdsn(gds_file_open, 'snp.id'))

# Samples to keep
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
            ifelse(is.null(args$include_samples),0,length(include_samples)), ifelse(is.null(args$include_samples),0,length(intersect(G_samples,include_samples))),
            ifelse(is.null(args$exclude_samples),0,length(exclude_samples)), ifelse(is.null(args$exclude_samples),0,length(intersect(G_samples,exclude_samples))),
            ifelse(is.null(args$include_samples)|is.null(args$exclude_samples),0,length(intersect(include_samples,exclude_samples))),
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

# Subset of matrix by samples
G_sub <- readex.gdsn(G, list(G_keep_samples, rep(TRUE, length(G_snps))))
# Transpose -> features x samples
G_sub <- t(G_sub)
# Set dim. names
dimnames(G_sub) <- list(G_snps,G_samples[G_keep_samples])
# Missing -> NA
G_sub[G_sub>2] <- NA # set missing values to NAs
if(args$verbose){ write('IMPORTANT: Missing value (genotype > 2, usually 3) was replaced ba NA',stdout()) }

# Save matrix
write.table(x=G_sub, file=args$bin_file, sep='\t', quote=FALSE, append=FALSE, row.names=TRUE, col.names=TRUE)
# Close GDS file
snpgdsClose(gds_file_open)
