#!/usr/bin/Rscript

RANKS <-  c('class',    'order',    'family',       'genus',    'species')
RANKS_ <- c('TaxClass','TaxOrder',  'TaxFamily',    'TaxGenus', 'TaxSpecies')

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))

## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--meta_rds', help='', required=TRUE, nargs='+')
    parser$add_argument('--obname', '-o', help='', required=TRUE)
    parser$add_argument('--ranks', default=RANKS, nargs='+', choices=RANKS)
    parser$add_argument('--src_path', help='')
    return(parser)
}

## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## Help functions
suppressMessages(source(sprintf('%s/utils.R',args$src_path)))

## Log file
log_file <- sprintf('%s.log', args$obname)
init_log(log_file, args)

## Data
tmp <- read_meta_rds(args$meta_rds)
meta <- tmp$meta
samples <- tmp$samples

## Tax. data
tax <- meta[,RANKS_] # select relevant columns
colnames(tax) <- RANKS # change column names
for(i in 1:ncol(tax)){ # convert taxa
    tax[,i] <- convert_taxa(tax[,i])
}

## Krona input data to create .xml
# Krona TEXT file
tax_text_file <- sprintf('%s.txt', args$obname)
# all unique species taxa
species_taxa <- sort(unique(tax$species))
# iterate over all taxa and save their counts and lineages
for(species_taxon in species_taxa){
    # get only samples assigned to that species
    tax_subset <-tax[tax$species == species_taxon, args$ranks]
    # lineage
    lineage <- tax_subset[!duplicated(tax_subset),,drop=FALSE] # the only unique row
    testit::assert(nrow(lineage)==1) # check
    # Abbr. genus name in species names if genus given
    if(!is.na(lineage[1,'genus'])){
        lineage[1,'species'] <- abbr_taxon(lineage[1,'species'])
    }
    # proc.
    lineage <- as.vector(na.omit(unlist(lineage[,args$ranks]))) # vector -> rm NAs -> vector

    write(sprintf('%d\t%s', nrow(tax_subset), paste(convert_taxa(lineage), collapse='\t')), file=tax_text_file, append=species_taxon != species_taxa[1])
}

## Call krona
# Krona HTML file
tax_html_file <- sprintf('%s.html', args$obname)
# CMD
cmd <- sprintf("ktImportText -n '' -o %s %s", tax_html_file, tax_text_file)

write(sprintf("CMD: %s", cmd), file=log_file, append=TRUE)
system(command=cmd)
