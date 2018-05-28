#!/usr/bin/Rscript

## NOTE IMPORTANT
##  Modified version of 10k_bacs_MIC_processing_v2.R
##      - Use given taxonomy lineages instead of calling taxize to find appropriate breakpoints
##      - Fix bug: Salmonella samples should use Salmonella breakpoints and not those of Enterobacteriaceae
## NOTE
##  Works with respect to KielID, i.e. all samples without KielID will be removed
##  Most things are hardcoded so be sure column names of relevant files are consistent with the code

## LIBS
suppressMessages(library(argparse)) # args parsing
suppressMessages(library(testit))   # assertions
suppressMessages(require(gtools))   # mixedsort

TAX_RANKS = c('Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum')

## FUNCTIONS
## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('-tax_file',    help='Should contain KielIDs and taxonomy lineage (species, genus, family, order, class, phylum)', required=TRUE)
    parser$add_argument('-tax_id_col',  help='Which column (name) in taxonomy file contains the sample KielID', required=TRUE)
    parser$add_argument('-id_file',     help='ID mapping file', required=TRUE)
    parser$add_argument('-mic_file',    help='Raw MIC file', required=TRUE)
    parser$add_argument('-mic_bp_file', help='MIC breakpoint file', required=TRUE)
    parser$add_argument('-o_dir',       help='Output directory', required=TRUE)
    parser$add_argument('-o_bname',     help='Output basename (without extension)', required=TRUE)
    
    parser$add_argument('--subset',     help='Create sub-matrix of binary res. profiles, e.g. \"Salmonella enterica:CP\"', nargs='+')
    
    parser$add_argument('--verbose', help='Print additional information', action='store_true')
    return(parser)
}

## Read ID mappig file
read_ids <- function(id_file, verbose=TRUE){
    ids <- read.csv(file=id_file, header=FALSE, sep='\t', stringsAsFactors=FALSE)
    colnames(ids) <- c('KielID','SacramentoID')
    if(verbose){ write(sprintf("\nRead in ID mapping file %s",id_file), stdout()); print(head(ids)) }
    return(ids)
}

## Read MIC data file
read_mics <- function(mic_file, verbose=TRUE){
    mics <- read.csv(file=mic_file, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
    # change one column name
    colnames(mics)[colnames(mics)=='Isolate'] <- 'SacramentoID'
    # only ID column and MICs
    mics <- mics[, grepl(' MIC$', colnames(mics)) | grepl('SacramentoID', colnames(mics))]
    # remove "MIC" from column names
    colnames(mics) <- sapply(colnames(mics), function(n){ ifelse(n=='SacramentoID',n,gsub(pattern='[[:blank:]]{1}MIC$',replacement='',x=n)) })
    
    if(verbose){ write(sprintf("\nRead in MIC data file %s",mic_file), stdout()); print(head(mics)) }
    return(mics)
}

## Read in MIC breakpoint file
read_bps <- function(bp_file, verbose=TRUE){
    bps <- read.csv(file=bp_file, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE, na.strings=c('NA','-','IE'))
    if(verbose){ write(sprintf("\nRead in MIC breakpoint file %s",bp_file), stdout()) }
    return(bps)
}

## Read in sample taxonomy file
read_tax <- function(tax_file, tax_id_col, tax_ranks=TAX_RANKS, verbose=TRUE){
    tax <- read.csv(file=tax_file, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    # check columns
    for(tax_rank in tax_ranks){
        testit::assert(sprintf('Rank %s not found in column names of taxonomy table',tax_rank), tax_rank %in% colnames(tax))
    }
    # change ID column name
    colnames(tax)[colnames(tax)==tax_id_col] <- 'KielID'
    # set row names
    rownames(tax) <- tax$KielID
    
    if(verbose){ write(sprintf("\nRead in taxonomy file %s",tax_file), stdout()); print(head(tax)) }
    return(tax)
}

## Create taxonomy lineages table
get_tax_lineages <- function(tax, tax_ranks=TAX_RANKS, verbose=TRUE){
    # Select only relevant columns
    tax_tl <- tax[,tax_ranks]
    # Keep only unique rows
    tax_tl <- unique(tax_tl)
    # Set row names
    rownames(tax_tl) <- tax_tl[,'Species']
    
    if(verbose){ write("\nTaxonomy lineage table:", stdout()); print(head(tax_tl)) }
    return(tax_tl)
}

## Process MIC value
## At the moment: Ignore the second component concentration 
process_mic <- function(mic){
    if(is.na(mic)){ return(NA) }
    # split MIC if two values
    mic <- unlist(strsplit(mic,'/'))
    testit::assert(length(mic)==1 | length(mic)==2)
    # remove sign, cast to numeric
    mic_ <- as.numeric( gsub(pattern='(<=|>)', replacement='', x=mic) )
    # process if sign
    if(grepl(pattern='<=',x=mic[1])){
        mic_ <- mic_/2
    }
    else if(grepl(pattern='>',x=mic[1])){
        mic_ <- mic_*2 ## incorrect for P/T (T conc. fixed at 4)
    }
    return(mic_[1]) # return only first/left value
}

## Process MIC table
process_mics <- function(mic_tab, ids, verbose=TRUE){
    # Sort columns by column names
    mics <- mics[,sort(colnames(mics))]
    # Map to KielIDs, NA if no match found
    mics$KielID <- sapply(mics$SacramentoID, function(sID){ ids$KielID[ match(x=sID, table=ids$SacramentoID) ] })
    # Remove samples without KielID
    if(verbose){
        write(sprintf(
            "From %d samples in MIC table %d have no KielID: %s",
            nrow(mics),sum(is.na(mics$KielID)),paste(mics$SacramentoID[is.na(mics$KielID)],collapse=';')
        ), stdout())
    }
    mics <- mics[!is.na(mics$KielID),]  
    # Keep only MIC columns
    rownames(mics) <- mics$KielID
    mics <- mics[,! colnames(mics) %in% c('KielID','SacramentoID') ]
    testit::assert('There are no 21 columns!',ncol(mics)==21)
    # Process
    mics <- apply(mics,c(1,2),function(mic){ process_mic(mic) })
    # Add samples with KielID not contained in MICs
    if(verbose){ write(sprintf("From %d samples in ID table %d have no MICs",nrow(ids),length(setdiff(ids$KielID,rownames(mics)))), stdout()) }
    for(kID in setdiff(ids$KielID,rownames(mics))){
        mics <- rbind(mics, NA); rownames(mics)[nrow(mics)] <- kID
    }
    # Sort by KielID
    mics <- mics[gtools::mixedsort(rownames(mics), decreasing=TRUE),]
    return(mics)
}

## Taxon mapping
get_taxon_mapping <- function(taxa, tax_ln, orgs, verbose=TRUE){
    if(verbose){ write('Matching taxa to taxa in MIC breakpoint file...', stdout()) }
    
    map <- data.frame( Taxon=taxa, bpTaxon=NA, row.names=taxa, stringsAsFactors=FALSE, check.names=FALSE )
    
    # map to MIC bp organism
    for(taxon in taxa){
        # taxon lineage
        taxon_tl <- tax_ln[taxon, , drop=FALSE]
        # from lowest to highest rank
        for(tax_rank in colnames(taxon_tl)){
            taxon_taxon <- unlist(taxon_tl[1, tax_rank])
            if(!is.na(taxon_taxon) & taxon_taxon %in% orgs){
                map[taxon,'bpTaxon'] <- taxon_taxon
                break
            }
        }
        if(verbose){ write(sprintf("\t%s: %s",taxon,map[taxon,'bpTaxon']), stdout()) }
    }
    return(map)
}

## Classify MICs
class_mics <- function(mics, tax, tax_map, bps, verbose=TRUE){
    mics_SIR <- matrix(NA, nrow=nrow(mics), ncol=ncol(mics), dimnames=dimnames(mics))
    bps_SIR  <- matrix(NA, nrow=nrow(mics), ncol=ncol(mics), dimnames=dimnames(mics))
    
    for(kID in rownames(mics)){
        taxon    <- tax[kID,'Species'] # sample taxon
        taxon_bp <- tax_map[taxon,'bpTaxon'] # mapped breakpoint organism
        if(is.na(taxon_bp)){
            write(sprintf("%s (%s): no mapped organism for MIC breakpoints", kID, taxon), stdout()) 
            next
        }
        for(drug in colnames(mics)){
            mic  <- mics[kID, drug] # sample MIC value
            R_bp <- bps[bps$Drug_abbr==drug,grepl(sprintf("%s_R",gsub(' ','_',taxon_bp)),colnames(bps))] # breakpoint for "resistant"
            S_bp <- bps[bps$Drug_abbr==drug,grepl(sprintf("%s_S",gsub(' ','_',taxon_bp)),colnames(bps))] # breakpoint for "susceptible"
            
            if(is.na(R_bp) | is.na(S_bp)){ ## No breakpoints
                testit::assert(is.na(R_bp) & is.na(S_bp))
            }
            else{
                bps_SIR[kID,drug] <- paste(S_bp,R_bp,sep=';')
                if(is.na(mic)){
                    mics_SIR[kID,drug] <- NA
                }
                else if(mic > R_bp){ ## EUCAST: >
                    mics_SIR[kID,drug] <- 'R'
                }
                else if(mic <= S_bp){ ## EUCAST: <=
                    mics_SIR[kID,drug] <- 'S'
                }
                else{
                    mics_SIR[kID,drug] <- 'I'
                }
            }
            if(verbose){
                write(sprintf(
                    "%s (%s): %s MIC = %s -> BP for %s for drug %s: S=%s, R=%s --> %s",
                    kID, taxon, drug, as.character(mic),
                    taxon_bp, drug, as.character(S_bp), as.character(R_bp), mics_SIR[kID,drug]
                ), stdout()) 
            }
        }
    }
    return(list(sir=mics_SIR,bps=bps_SIR))
}



## ARGS
# test_args <- c('-tax_file ~/Bacteria_project/10k_bacs/Sample_data//kraken_tax.csv -tax_col Species -o_dir ~/Bacteria_project/10k_bacs --verbose')
# args <- get_argparser()$parse_args(strsplit(test_args,' ')[[1]])
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

# Print info
if(args$verbose){
    write(sprintf(
        "%s\n%s\n\nArguments:\n%s",
        format(Sys.time(), "%Y.%m.%d, %H:%M:%S"),
        paste(commandArgs(),collapse=' '),
        paste(paste('\t',names(unlist(args)),unlist(args),sep=' '),collapse='\n')
    ), stdout())
}

## DATA
ids <- read_ids(id_file=args$id_file, verbose=args$verbose)
mics<- read_mics(mic_file=args$mic_file, verbose=args$verbose)
bps <- read_bps(bp_file=args$mic_bp_file, verbose=args$verbose)
tax <- read_tax(tax_file=args$tax_file, tax_id_col=args$tax_id_col, tax_ranks=TAX_RANKS, verbose=args$verbose)

# Taxonomy lineages
tax_ln <- get_tax_lineages(tax=tax, tax_ranks=TAX_RANKS, verbose=args$verbose)

# Process MICs
mics2 <- process_mics(mic_tab=mics, ids=ids, verbose=args$verbose)

# Find for each sample the corresponding organism for the breakpoints
tax_map <- get_taxon_mapping(
    taxa=sort(unique(tax$Species)),
    tax_ln=tax_ln,
    orgs=sort(unique( sapply(colnames(bps),function(x){ ifelse(grepl('_(S|R)$',x),gsub('_',' ',gsub('_(S|R)$','',x)),NA) }) )),
    verbose=args$verbose
)

# Classify as S/I/R
SIR <- class_mics(mics=mics2, tax=tax, tax_map=tax_map, bps=bps, verbose=TRUE)
BP  <- SIR[['bps']]
SIR <- SIR[['sir']]

## SAVE
write.table(x=mics2, file=sprintf("%s/%s_prMICs.csv",args$o_dir,args$o_bname), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=BP,    file=sprintf("%s/%s_SR_BP.csv",args$o_dir,args$o_bname),  sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=SIR,   file=sprintf("%s/%s_SIR.csv",args$o_dir,args$o_bname),    sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(
    x=apply(SIR, c(1,2), function(x){ ifelse(is.na(x),NA,as.numeric(x=='R')) }),
    file=sprintf("%s/%s_R.csv",args$o_dir,args$o_bname),
    sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE
)

if(!is.null(args$subset)){
    for(s in args$subset){
        s_ = unlist(strsplit(s,':'))
        s_sp <- s_[1]
        s_dr <- ifelse(length(s_)==1,colnames(SIR),s_[2])

        testit::assert(sprintf('Species %s not in taxonomy data', s_sp), s_sp %in% tax$Species)
        for(d in s_dr){
            testit::assert(sprintf('Drug %s not in created resistance table', d), d %in% colnames(SIR))
        }

        if(args$verbose){
            write(sprintf("Creating sub-matrix of binary res. profiles for species %s and drugs %s", s_sp, paste(s_dr, collapse=', ')), stdout())
        }
        write.table(
            x=apply(SIR[rownames(tax)[tax$Species==s_sp],s_dr,drop=FALSE], c(1,2), function(x){ ifelse(is.na(x),NA,as.numeric(x=='R')) }),
            file=sprintf("%s/%s_%s_%s_R.csv",args$o_dir,args$o_bname, gsub(' ','_',s_sp), ifelse(length(s_dr)==1,s_dr,'all')),
            sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE
        )
    }
}
