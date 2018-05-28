#!/usr/bin/Rscript

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(tools))
suppressMessages(library(taxize))


## ARGS
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('-tax_file',help='Processed and merged MALDI results', required=TRUE)
    parser$add_argument('-tax_col',help='', required=TRUE)
    parser$add_argument('-id_col',help='', required=TRUE)
    parser$add_argument('-o_file_tax',help='Output file for full taxonomy table', required=TRUE)
    parser$add_argument('-o_bname_krona',help='Path/basename for the Krona files', required=TRUE)
    parser$add_argument('-krona_script',help='', default='Utils/create_krona_xml.R')
    parser$add_argument('--verbose', help='Print additional information', action='store_true')
    
    return(parser)
}

## Read in sample taxonomy file
read_tax <- function(tax_file, tax_col, id_col, verbose=TRUE){
    tax <- read.csv(file=tax_file, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    
    tax <- tax[,c(id_col,tax_col)]
    colnames(tax) <- c('KielID','Taxon')
    rownames(tax) <- tax$KielID
    
    if(verbose){ write(sprintf("\nRead in taxonomy file %s",tax_file), stdout()); print(head(tax)) }
    return(tax)
}

## FILES AND DATA
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## DATA
tax <- read_tax(tax_file=args$tax_file, tax_col=args$tax_col, id_col=args$id_col, verbose=args$verbose)

## Processing
tax_table <- data.frame(
    KielID=tax$KielID,
    Organism=tax$Taxon,
    Species=NA, Genus=NA, Family=NA, Order=NA, Class=NA, Phylum=NA,
    stringsAsFactors=FALSE
)

orgs <- sort(unique(tax_table$Organism))
source('Utils/taxize_with_DB.R')
tax_map <- get_tax_simple(orgs=orgs, db='ncbi', use_pb=TRUE)

for(i in 1:nrow(tax_table)){
    orgs <- tax_table$Organism[i]
    
    if(is.na(orgs)){ next }
    
    orgs <- strsplit(orgs, ';')[[1]] # if more than one organism (separated by ";")
    org_species <- org_genus <- org_family <- org_order <- org_class <- org_phylum <- NULL
    
    for(org in orgs){
        if(!is.na(tax_map[org,'Species'])){ org_species <- c(org_species,tax_map[org,'Species']) }
        org_genus  <- c(org_genus, tax_map[org,'Genus'])
        org_family <- c(org_family,tax_map[org,'Family'])
        org_order  <- c(org_order, tax_map[org,'Order'])
        org_class  <- c(org_class, tax_map[org,'Class'])
        org_phylum <- c(org_phylum,tax_map[org,'Phylum'])
    }
    
    # Species
    tax_table$Species[i] <- paste(unique(org_species),collapse='/')
    tax_table$Genus[i]   <- paste(unique(org_genus),collapse='/')
    tax_table$Family[i]  <- paste(unique(org_family),collapse='/')
    tax_table$Order[i]   <- paste(unique(org_order),collapse='/')
    tax_table$Class[i]   <- paste(unique(org_class),collapse='/')
    tax_table$Phylum[i]  <- paste(unique(org_phylum),collapse='/')
}

# Save
tax_file <- args$o_file_tax
write.table(x=tax_table, file=tax_file, row.names=FALSE, sep='\t', quote=FALSE)

## Create xml
xml_obname  <-  args$o_bname_krona
ranks    <- c('Species','Genus','Family','Order','Class')
cmd <- sprintf('Rscript %s %s %s --ranks %s --verbose',args$krona_script,tax_file,xml_obname,paste(ranks,collapse=' '))
if(args$verbose){ write(cmd, stdout()) }
system(command=cmd)

## Call krona
cmd <- sprintf("ktImportXML -o %s %s",sprintf("%s.html",xml_obname),sprintf("%s.xml",xml_obname))
if(args$verbose){ write(cmd, stdout()) }
system(command=cmd)
