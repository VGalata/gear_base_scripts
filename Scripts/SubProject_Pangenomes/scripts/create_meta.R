#!/usr/bin/Rscript

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))

## Help functions
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--db_tab', help='db_samples_*', required=TRUE)
    parser$add_argument('--smeta', help='saureus_meta.tsv', required=TRUE)
    parser$add_argument('--usmap', help='us_max_locs_map.tsv', required=TRUE)
    parser$add_argument('--samples', help='all_samples.txt', required=TRUE)
    parser$add_argument('--obname', required=TRUE)
    parser$add_argument('--src_path', required=TRUE)
    return(parser)
}

# state data
state_data <- data.frame(
    Name=state.name,
    Abbr=state.abb,
    Region=as.character(state.region),
    row.names=state.name, check.names=FALSE, stringsAsFactors=FALSE
)

usstate2region <- function(x){ return(state_data[x,'Region']) }

usstate2abbr   <- function(x){ return(state_data[x,'Abbr']) }

## add US additional mapping to locations
read_us_states_mapping <- function(map_file='/home/vgalata/git_repos/Bacteria/Data/pangenomes/us_max_locs_map.tsv'){
    return(
        mapping <- read.csv(file=map_file, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    )
}

add_us_states <- function(df, mapping, loc_col='Max_Name', min_loc_col='Min_Name'){
    # checks
    testit::assert(all(sapply(c('Name', 'ShortName', 'LongName'), function(x){ x %in% colnames(mapping) })))
    # map
    for(i in 1:nrow(df)){
        loc_name <- df[i, loc_col]
        if(loc_name %in% mapping$Name){
            # check
            testit::assert(df[i,min_loc_col] == 'US')
            # match
            j <- which(mapping$Name == loc_name)
            # long/short name
            df$USStateShort[i] <- mapping$ShortName[j]
            df$USStateLong[i]  <- mapping$LongName[j]
            df$USStateRegion[i]<- usstate2region(df$USStateLong[i])
        }
    }
    return(df)
}

## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## Help functions
suppressMessages(source(sprintf('%s/utils.R', args$src_path)))

## Log
log_file <- sprintf('%s.log', args$obname)
init_log(log_file, args)

## Data
tab1 <- read.csv(file=args$db_tab, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
tab2 <- read.csv(file=args$smeta, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
usmap<- read_us_states_mapping(args$usmap)
ids  <- unlist(read.csv(file=args$samples, header=FALSE, stringsAsFactors=FALSE))

write(sprintf(
    'DB data table %d x %d\nS. aureus meta data table %d x %d\nUS states map %d x %d\nIDs %d',
    nrow(tab1), ncol(tab1),
    nrow(tab2), ncol(tab2),
    nrow(usmap), ncol(usmap),
    length(ids)
    ),
    log_file, append=TRUE
)

## Checks
# only unique IDs
testit::assert(length(ids) == length(unique(ids)))
# all IDs in 1st table
testit::assert(all(sapply(ids, function(x){ x %in% tab1$sample___seq_id })))
# complete intersection between two tables
testit::assert(length(intersect(tab1$sample___seq_id, tab2$Kiel_NGS_ID)) == nrow(tab2))

## Select columns in 1st table
df <- data.frame(
    # ID (Kiel NGS ID)
    ID=tab1$sample___seq_id,
    # Which dataset TODO
    DataSet=NA,
    # Collection/isolation date
    Date=as.Date(tab1$sample___date, '%Y.%m.%d'),
    # Source (where collected for 10k bac.s, study source for S. aureus)
    Source=tab1$sample___source__name,
    # Location (min = country/continent, max = e.g. state)
    MinLoc=tab1$sample___source__loc_min__name,
    MaxLoc=tab1$sample___source__loc_max__name,
    # US state assignment (short and long name, region) TODO
    USStateShort=NA, USStateLong=NA, USStateRegion=NA,
    # Taxonomy lineage
    TaxPhylum=tab1$taxlin___phylum,
    TaxClass=tab1$taxlin___class,
    TaxOrder=tab1$taxlin___order,
    TaxFamily=tab1$taxlin___family,
    TaxGenus=tab1$taxlin___genus,
    TaxSpecies=tab1$taxlin___species,
    # Taxonomy quality
    TaxSens=tab1$sample___tax_assign_sample__pct_sens,
    TaxPrec=tab1$sample___tax_assign_sample__pct_prec,
    TaxUncl=tab1$sample___tax_assign_sample__pct_uncl,
    # Assembly quality
    AssemblyNum=tab1$sample___assembly_sample__aq_num,
    AssemblyLen=tab1$sample___assembly_sample__aq_len,
    AssemblyMaxLen=tab1$sample___assembly_sample__aq_max_len,
    AssemblyGC=tab1$sample___assembly_sample__aq_GC,
    AssemblyNs=tab1$sample___assembly_sample__aq_Ns,
    AssemblyN50=tab1$sample___assembly_sample__aq_N50,
    AssemblyL50=tab1$sample___assembly_sample__aq_L50,
    #
    row.names=tab1$sample___seq_id, stringsAsFactors=FALSE, check.names=FALSE
)

## Add dataset assignment
df[tab2$Kiel_NGS_ID,  'DataSet'] <- 'saureus'
df[is.na(df$DataSet), 'DataSet'] <- '10kbacs'

testit::assert(sum(df$DataSet=='10kbacs')==10086)
testit::assert(sum(df$DataSet=='saureus')==1001)

## Add S. aureus meta data
testit::assert(all(is.na(df[tab2$Kiel_NGS_ID,  'Date'])))
testit::assert(all(is.na(df[tab2$Kiel_NGS_ID,  'MinLoc'])))

df[tab2$Kiel_NGS_ID, 'Date']   <- as.Date(as.character(tab2$IsolationYear), '%Y')
df[tab2$Kiel_NGS_ID, 'MinLoc'] <- tab2$Country
df[tab2$Kiel_NGS_ID, 'Source'] <- tab2$StudySource

## Add US states info
df <- add_us_states(df=df, mapping=usmap, loc_col='MaxLoc', min_loc_col='MinLoc')

## De/Select samples
df_1 <- df[ids,]
df_2 <- df[setdiff(rownames(df),ids),]

testit::assert(all(sapply(ids, function(x){   x %in% rownames(df_1) })))
testit::assert(all(sapply(ids, function(x){ ! x %in% rownames(df_2) })))
testit::assert(length(intersect(rownames(df_1), rownames(df_2)))==0)
testit::assert(nrow(df_1)+nrow(df_2) == nrow(df))

## Save
write.table(x=df_1,  file=sprintf('%s.tsv', args$obname), sep='\t', row.names=FALSE, quote=FALSE)
saveRDS(object=df_1, file=sprintf('%s.rds', args$obname))

write.table(x=df_2,  file=sprintf('%s_rest.tsv', args$obname), sep='\t', row.names=FALSE, quote=FALSE)
saveRDS(object=df_2, file=sprintf('%s_rest.rds', args$obname))
