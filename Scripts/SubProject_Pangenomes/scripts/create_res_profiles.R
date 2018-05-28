#!/usr/bin/Rscript

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))

## Help functions
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--db_tab', help='db_samples_*', required=TRUE)
    parser$add_argument('--obname', required=TRUE)
    parser$add_argument('--src_path', required=TRUE)
    return(parser)
}

## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## Help functions
suppressMessages(source(sprintf('%s/utils.R', args$src_path)))

## Log
log_file <- sprintf('%s.log', args$obname)
init_log(log_file, args)

## Data
tab <- read.csv(file=args$db_tab, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
write(sprintf('DB data table %d x %d',nrow(tab), ncol(tab)), log_file, append=TRUE)

## Select columns in 1st table
df_mic <- data.frame(
    tab[,grepl('^drugres___.*___mic$',colnames(tab))],
    row.names=tab$sample___seq_id, stringsAsFactors=FALSE, check.names=FALSE
)
colnames(df_mic) <- gsub('drugres___', '', gsub('___mic', '', colnames(df_mic)))
write(sprintf('MIC data table %d x %d',nrow(df_mic), ncol(df_mic)), log_file, append=TRUE)

df_res <- data.frame(
    tab[,grepl('^drugres___.*___value$',colnames(tab))],
    row.names=tab$sample___seq_id, stringsAsFactors=FALSE, check.names=FALSE
)
colnames(df_res) <- gsub('drugres___', '', gsub('___value', '', colnames(df_res)))
write(sprintf('Res data table %d x %d',nrow(df_res), ncol(df_res)), log_file, append=TRUE)

## Save
write.table(x=df_mic,  file=sprintf('%s_mic.tsv', args$obname), sep='\t', quote=FALSE)
saveRDS(object=df_mic, file=sprintf('%s_mic.rds', args$obname))

write.table(x=df_res,  file=sprintf('%s_res.tsv', args$obname), sep='\t', quote=FALSE)
saveRDS(object=df_res, file=sprintf('%s_res.rds', args$obname))
