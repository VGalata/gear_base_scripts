#!/usr/bin/Rscript

save_session_info <- function(filename, append=TRUE){
    # write session info. to file
    write("Session info:", file=filename, append=append)
    write(capture.output(sessionInfo()), file=filename, append=append)
    write("\n", file=filename, append=append)
}

save_args <- function(args, filename, append=TRUE){
    write("Arguments:", file=filename, append=TRUE)
    for(a_name in sort(names(args))){
        write(sprintf('%s: %s', a_name, as.character(args[a_name])), file=filename, append=TRUE)
    }
    write("\n", file=filename, append=TRUE)
}

init_log <- function(filename, args){
    # time stamp
    write(capture.output(Sys.time()), file=filename, append=FALSE)
    # session information
    save_session_info(log_file, append=TRUE)
    # given arguments
    save_args(args, log_file, append=TRUE)
}

read_meta_rds <- function(meta_files){
    meta <- NULL
    samples <- list()
    for(meta_file in meta_files){
        tmp  <- readRDS(meta_file)
        meta <- rbind(meta, tmp)
        if(grepl('rest', meta_file)){
            samples$rest <- rownames(tmp)
        } else {
            samples$main <- rownames(tmp)
        }
    }
    return(list(meta=meta, samples=samples))
}

read_in_pam <- function(pam_file){
    return(read.csv(file=pam_file, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
}

taxon_from_filename <- function(filename){
    # basename
    taxon <- basename(dirname(filename))
    # rm extension
    # taxon <- unlist(strsplit(taxon, '\\.'))[1]
    # replace underscores
    taxon <- gsub('_', ' ', taxon)
    return(taxon)
}

read_pg_stats <- function(filename){
    return(read.csv(file=filename, header=TRUE, sep='\t', stringsAsFactors=FALSE, row.names=1))
}

abbr_taxon <- function(taxon){
    if(taxon=='Other'){ return(taxon) }
    taxon <- unlist(strsplit(taxon, ' '))
    return(sprintf('%s. %s', substr(taxon[1], 1, 1), paste(taxon[-1], collapse=' ')))
}

abbr_taxon2 <- function(taxon){
    if(taxon=='Other'){ return(taxon) }
    taxon <- unlist(strsplit(taxon, ' '))
    return(paste(sapply(taxon, function(x){ sprintf('%s.', substr(x, 1, 2)) }), collapse=' '))
}

abbr_taxa <- function(taxa, method=1){
    if(method == 1){
        return( sapply(as.character(taxa), abbr_taxon) )
    } else if(method == 2){
        return( sapply(as.character(taxa), abbr_taxon2) )
    } else {
        return( sapply(as.character(taxa), abbr_taxon) )
    }
}

# "<taxon name>:<taxon ID>" -> "<taxon name>"
convert_taxon <- function(taxon){
    taxon <- as.character(taxon)
    return(ifelse(is.na(taxon), NA, sub(':[0-9]+$', '', taxon)))
}

# remove taxonomy ID from taxon string "<name>:<id>"
convert_taxa <- function(taxa){
    return(
        sapply(as.character(taxa), convert_taxon)
    )
}

#
hmm_tblout_cols <- c(
    'target_name', 'target_accession',
    'query_name', 'query_accession',
    'full_sequence_Evalue', 'full_sequence_score', 'full_sequence_bias',
    'best_1_domain_Evalue', 'best_1_domain_score', 'best_1_domain_bias',
    'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target'
)

geneID2sampleID <- function(gID){
    return( unlist(strsplit(gID, '_'))[1] )
}

getYearFromDate <- function(d){
    return( as.numeric(format(d,'%Y')) )
}
