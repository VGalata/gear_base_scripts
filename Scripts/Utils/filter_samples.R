#!/usr/bin/Rscript

## NOTE
## Script for filtering samples based on:
##  1) Number of reads in raw/trimmed fastq.gz files
##  2) Assembly quality (QUAST output)
##  3) Taxnonomy analysis (Kraken results)
## Creates a summary table and sample lists based on taxnomic grouping



## LIBS
suppressMessages(library(argparse)) # args parsing
suppressMessages(library(testit))   # assertions
suppressMessages(library(taxize))   # taxonomy
# suppressMessages(library(tools))    # for string operations


## ARGS
# Command line argument parser
get_argparser <- function(){
    script_descr="Script for sample filtering based on assembly quality (default: NCBI RefSeq min. criteria) and taxonomy analysis"
    parser <- ArgumentParser(description=script_descr)
    
    # required parameters: input
    parser$add_argument('-a_file', help='Assembly quality: Required columns: Assembly, # contigs, N50, L50', required=TRUE)
    parser$add_argument('-k_file', help='Taxonomy summary: Required columns: ID, <rank>, <rank>_count, <rank>_rank_count, Unclassified_count and TotalCount', required=TRUE)
    # required parameters: output
    parser$add_argument('-o_dir',  help='Output directory', required=TRUE)
    
    # filtering criteria: assembly quality
    parser$add_argument('--max_contigs', help='Max. number of contigs, default is 1000 (NCBI RefSeq rule)',                             type="integer", default=1000)
    parser$add_argument('--min_N50',     help='Min. N50 value, default is 5000 (NCBI RefSeq rule)',                                     type="integer", default=5000)
    parser$add_argument('--max_L50',     help='Max. L50 value, default is 200 (NCBI RefSeq rule)',                                      type="integer", default=200)
    
    # filtering criteria: taxnomy
    parser$add_argument('--rank',        help='Which rank was used to create the Kraken summary file, default is Species',              default="Species")
    parser$add_argument('--min_sens',    help='Minimum (>=) sensitivity=(# mapped to taxon)/(# total), default is 50',                  type="integer", default=50)
    parser$add_argument('--min_prec',    help='Minimum (>=) precision=(# mapped to taxon)/(# mapped at taxon level), default is 75',    type="integer", default=75)
    parser$add_argument('--max_uncl',    help='Maximum (<=) percentage of unclassified , default is 30',                                type="integer", default=30)
    parser$add_argument('--shigella_min_sens', help='Minimum (>=) sensitivity for Shigella spp., default is 0',                         type="integer", default=0)
    parser$add_argument('--shigella_min_prec', help='Minimum (>=) precision  for Shigella spp., default is 60',                         type="integer", default=60)
    parser$add_argument('--shigella_max_uncl', help='Maximum (<=) percentage for Shigella spp., default is 30',                         type="integer", default=30)
    
    # rest
    parser$add_argument('--min_num', help='Minimal number of samples (after filtering) to create sample lists, default is 50',          type="integer", default=50)
    parser$add_argument('--add_tax', help='Add further taxonomic information as genus, family etc. (using taxize, NCBI)',               action='store_true')
    
    # other
    parser$add_argument('--verbose', help='Print additional information', action='store_true')
    
    return(parser)
}

# Parse args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))
log_file <- sprintf("%s/kraken_sum_tax.log",args$o_dir)

# Print info
if(args$verbose){
    write(sprintf('\nGiven arguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
    write(sprintf('\n'),stdout())
}
# Log
write(sprintf(
    "%s\n%s\nArguments:\n%s",
    format(Sys.time(), "%Y.%m.%d, %H:%M:%S"),
    paste(commandArgs(),collapse=' '),
    paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n')
), log_file, append=FALSE)


## DATA

## Assembly quality
assem_quali <- read.csv(file=args$a_file, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
rownames(assem_quali) <- assem_quali$Assembly

if(args$verbose){ write(sprintf("> Assembly quality: %d samples",nrow(assem_quali)), stdout()) }
write(sprintf("> Assembly quality: %d samples",nrow(assem_quali)), log_file, append=TRUE)

## Taxonomy
taxonomy <- read.csv(file=args$k_file, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
rownames(taxonomy) <- taxonomy$ID

if(args$verbose){ write(sprintf("> Taxonomy: %d samples",nrow(taxonomy)), stdout()) }
write(sprintf("> Taxonomy: %d samples",nrow(taxonomy)), log_file, append=TRUE)

# Sensitivity and precision
taxonomy$Sens <- 100 * taxonomy[,sprintf("%s_count",args$rank)] / taxonomy$TotalCount
taxonomy$Prec <- 100 * taxonomy[,sprintf("%s_count",args$rank)] / taxonomy[,sprintf("%s_rank_count",args$rank)]

# Unclassified percentage
taxonomy$UnclassPct <- 100 * taxonomy$Unclassified_count / taxonomy$TotalCount

# Checks
if(sprintf('%s_pct',args$rank) %in% colnames(taxonomy)){
    testit::assert(all( abs(taxonomy$Species_pct - taxonomy$Sens) < 1e-1 ))
}
if('Unclassified_pct' %in% colnames(taxonomy)){
    testit::assert(all( abs(taxonomy$Unclassified_pct - taxonomy$UnclassPct) < 1e-1 ))
}

## Samples
samples <- sort( intersect( rownames(taxonomy),rownames(assem_quali) ) )

log_text <- sprintf("> There are %d samples in all input files...",length(samples))
if(args$verbose){ write(log_text, stdout()) }
write(log_text, log_file, append=TRUE)

assem_quali <- assem_quali[samples,]
taxonomy    <- taxonomy[samples,]


## FILTER

## By assembly quality
keep_by_aq <- sapply(1:nrow(assem_quali), function(i){
    if(any(is.na(assem_quali[i,]))){ return(FALSE) }
    return(
        assem_quali[i,'# contigs'] <= args$max_contigs &
        assem_quali[i,'N50'] >= args$min_N50 &
        assem_quali[i,'L50'] <= args$max_L50
    )
})
keep_by_aq <- samples[keep_by_aq]

log_text <- sprintf("> Filter by assembly quality: %d from %d were kept",length(keep_by_aq),length(samples))
write(log_text, log_file, append=TRUE)
if(args$verbose){ write(log_text, stdout()) }

## By Taxonomy
# Filter all samples
keep_by_ta <- sapply(1:nrow(taxonomy), function(i){
    taxon<- taxonomy[i, args$rank]    # Taxon
    sens <- taxonomy[i, 'Sens']       # Sensitivity
    prec <- taxonomy[i, 'Prec']       # Precision
    uncl <- taxonomy[i, 'UnclassPct'] # Unclassified pct.
    if(grepl('Shigella',taxon)){ # Shigella
        return( sens>=args$shigella_min_sens & prec>=args$shigella_min_prec & uncl<=args$shigella_max_uncl )
    }
    else{ # Not Shigella
        return( sens>=args$min_sens          & prec>=args$min_prec          & uncl<=args$max_uncl )
    }
})
keep_by_ta <- samples[keep_by_ta]

log_text <- sprintf("> Filter by taxonomy: %d from %d were kept",length(keep_by_ta),length(samples))
write(log_text, log_file, append=TRUE)
if(args$verbose){ write(log_text, stdout()) }

## Intersection
keep <- sort( intersect( keep_by_aq , keep_by_ta ) )

log_text <- sprintf("> Filter by all: %d from %d were kept",length(keep),length(samples))
write(log_text, log_file, append=TRUE)
if(args$verbose){ write(log_text, stdout()) }


## ADD TAX
tax_map <- NULL
if(args$add_tax){
    if(args$verbose){ write("> Retrieving taxonomy... This may take some time...", stdout()) }
    
    # All taxons appearing in Kraken table
    orgs <- sort(unique(taxonomy[,args$rank]))
    # Get additional taxonomy information
    source('Utils/taxize_with_DB.R')
    tax_map <- get_tax_simple(orgs=orgs, db='ncbi', use_pb=args$verbose)
    
    if(args$verbose){ write("\t> Done.", stdout()) }
}


## OUTPUT
result <- data.frame(
    KielID          =samples, # KielID
    Passed          =sapply(samples, function(kid){ kid %in% keep }), # If passed th filter or not
    Taxon           =taxonomy[samples,args$rank], # Taxon
    TaxonCount      =taxonomy[,sprintf("%s_count",args$rank)],
    TaxonRankCout   =taxonomy[,sprintf("%s_rank_count",args$rank)],
    UnclassCount    =taxonomy$Unclassified_count,
    TotalCount      =taxonomy$TotalCount,
    Sens            =round(taxonomy$Sens,2),
    Prec            =round(taxonomy$Prec,2),
    UnclassPct      =round(taxonomy$UnclassPct,2),
    Contigs         =assem_quali[samples,'# contigs'],
    N50             =assem_quali[samples,'N50'],
    L50             =assem_quali[samples,'L50'],
    stringsAsFactors=FALSE, check.names=FALSE, row.names=samples
)

# add additional tax. if asked
if(args$add_tax){
    result <- data.frame(
        result,
        TaxSpecies=sapply(result$Taxon, function(taxon){ tax_map[tax_map$Organism==taxon,'Species'] }),
        TaxGenus=  sapply(result$Taxon, function(taxon){ tax_map[tax_map$Organism==taxon,'Genus'] }),
        TaxFamily= sapply(result$Taxon, function(taxon){ tax_map[tax_map$Organism==taxon,'Family'] }),
        TaxOrder=  sapply(result$Taxon, function(taxon){ tax_map[tax_map$Organism==taxon,'Order'] }),
        TaxClass=  sapply(result$Taxon, function(taxon){ tax_map[tax_map$Organism==taxon,'Class'] }),
        TaxPhylum= sapply(result$Taxon, function(taxon){ tax_map[tax_map$Organism==taxon,'Phylum'] }),
        stringsAsFactors=FALSE, check.names=FALSE, row.names=rownames(result)
    )
}


## SAVE
## Tables
write.table(x=result, file=sprintf("%s/samples_all.csv", args$o_dir), sep='\t', row.names=FALSE, quote=FALSE)
write.table(x=result[result$Passed,], file=sprintf("%s/samples.csv", args$o_dir), sep='\t', row.names=FALSE, quote=FALSE)

## Filtered samples
write(x=setdiff(samples, keep_by_aq), file=sprintf("%s/samples_removed_aq.txt", args$o_dir), ncolumns=1, sep='\t', append=FALSE)
write(x=setdiff(samples, keep_by_ta), file=sprintf("%s/samples_removed_ta.txt", args$o_dir), ncolumns=1, sep='\t', append=FALSE)

## Sample lists
if(args$verbose){ write("> Tables were saved\n> Creating sample lists...", stdout()) }

# all samples: unfiltered, filtered
write(x=result$KielID, file=sprintf("%s/samples_all.txt", args$o_dir), ncolumns=1, sep='\t', append=FALSE)
write(x=result$KielID[result$Passed], file=sprintf("%s/samples.txt", args$o_dir), ncolumns=1, sep='\t', append=FALSE)

# grouping
n <- max(sapply(unique(result$Taxon), nchar)) # for printing
for(taxon in sort(unique(result$Taxon))){
    if(is.na(taxon)){ next } # should not happen
    
    # All samples of that taxon not filtered out -> KielID
    kIDs <- result$KielID[result$Passed & result$Taxon==taxon]
    
    # Print info
    if(args$verbose){
        m <- nchar(taxon); empty_str <- paste(rep(' ',1+n-m),collapse='') # for printing
        s <- sprintf( "\t%s%s--> %s: %d", taxon, empty_str, ifelse(length(kIDs)>=args$min_num,'create list','no list (too few samples)'), length(kIDs) )
        write(s,stdout())
    }
    write(s,log_file,append=TRUE)
    
    # Create list if enough samples
    if(length(kIDs)>=args$min_num){
        write(x=kIDs, file=sprintf("%s/samples_%s.txt", args$o_dir, gsub(pattern=' ',replacement='_',x=taxon)), ncolumns=1, sep='\t', append=FALSE)
    }
}
