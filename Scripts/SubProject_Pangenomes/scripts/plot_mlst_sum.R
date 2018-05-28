#!/usr/bin/Rscript

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))

# args <- list(
#     meta_rds='sample_data.rds',
#     mlst_tab='/share/runs/assembly/11k_summaries/mlst_tseemann_10kbacs_saureus_2.tsv',
#     src_path="~/git_repos/Bacteria/SubProject_Pangenomes/scripts"
# )


## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--meta_rds', help='', required=TRUE)
    parser$add_argument('--mlst_tab', help='', required=TRUE)
    parser$add_argument('--pdf', help='', required=TRUE)
    parser$add_argument('--src_path', help='')
    parser$add_argument('--topN', help='', type="integer", default=5)
    return(parser)
}

## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## Help functions
suppressMessages(source(sprintf('%s/utils.R',args$src_path)))

## Log file
log_file <- sprintf('%s.log', args$pdf)
init_log(log_file, args)

## Data
tmp     <- read_meta_rds(args$meta_rds)
meta    <- tmp$meta
samples <- tmp$samples$main
mlst    <- read.csv(file=args$mlst_tab, sep='\t', header=TRUE, stringsAsFactors=FALSE)
# check IDs
testit::assert(all(samples %in% mlst$KielID))
# process
rownames(mlst) <- mlst$KielID
mlst <- mlst[samples,]
# NAs
mlst$MLSTst[is.na(mlst$MLSTst)] <- 'unknown'
mlst$MLSTscheme[is.na(mlst$MLSTscheme)] <- 'unknown'
# check taxa
testit::assert(all(mlst$KrakenTaxon == convert_taxa(meta[samples,'TaxSpecies'])))

##
mlst_aggr <- aggregate(
    1:nrow(mlst),
    by=list(Taxon=mlst$KrakenTaxon, ST=mlst$MLSTst, Scheme=mlst$MLSTscheme),
    FUN=function(x){ length(x) }, drop=TRUE
)

## Plots
pdf(args$pdf, width=6, height=6)
for(taxon in sort(unique(mlst_aggr$Taxon))){
    # filter
    taxon_mlst <- mlst_aggr[mlst_aggr$Taxon==taxon,]
    # sort by count (no unknowns)
    taxon_mlst_ <- taxon_mlst[taxon_mlst$ST != "unknown" & taxon_mlst$Scheme != "unknown",,drop=FALSE]
    taxon_mlst_ <- taxon_mlst_[with(taxon_mlst_, order(-x)), ] # big -> small
    # discarded
    if(nrow(taxon_mlst_) > args$topN){
        discarded <- rownames(taxon_mlst_)[(args$topN+1):nrow(taxon_mlst_)]
        taxon_mlst_[discarded, 'ST'] <- ''
        taxon_mlst_[discarded, 'Scheme'] <- 'other'
    }
    # top N + discarded + unknowns
    taxon_mlst <- rbind(
        taxon_mlst_,
        taxon_mlst[taxon_mlst$ST == "unknown" | taxon_mlst$Scheme == "unknown",]
    )
    # aggregate discarded
    taxon_mlst <- aggregate(taxon_mlst$x, by=list(ST=taxon_mlst$ST, Scheme=taxon_mlst$Scheme), FUN=sum)
    # plot labels
    taxon_mlst$Label <- paste(taxon_mlst$Scheme, taxon_mlst$ST, sep='/')
    taxon_mlst$Label <- gsub('/$', '', taxon_mlst$Label) # if discarded
    taxon_mlst$Label <- gsub('unknown/unknown$', 'unknown', taxon_mlst$Label) # if no scheme/ST
    taxon_mlst$Label <- gsub('/unknown$', '/new', taxon_mlst$Label) # <scheme>/unknown -> <scheme>/new

    # check the number of samples
    testit::assert( sum(taxon_mlst$x) == sum(mlst$KrakenTaxon == taxon) )
    # log
    write(sprintf('\n%s: %d', taxon, sum(taxon_mlst$x)), log_file, append=TRUE)
    write(sprintf('Proportion of "other": %.1f', 100 * taxon_mlst$x[taxon_mlst$Label == "other"] / sum(taxon_mlst$x)), log_file, append=TRUE)
    write(sprintf('Proportion of "unknown": %.1f', 100 * taxon_mlst$x[taxon_mlst$Label == "unknown"] / sum(taxon_mlst$x)), log_file, append=TRUE)
    write(capture.output(print(taxon_mlst)), log_file, append=TRUE)

    scale_labels <- setdiff(taxon_mlst$Label, c('other', 'unknown'))
    scale_values <- colorRampPalette(c("#CCFF99", "#00CC66"))(length(scale_labels))
    if('other' %in% taxon_mlst$Label){
        scale_labels <- c(scale_labels, 'other')
        scale_values <- c(scale_values, "#3399FF")
    }
    if('unknown' %in% taxon_mlst$Label){
        scale_labels <- c(scale_labels, 'unknown')
        scale_values <- c(scale_values, "#FFCCCC")
    }

    p <-
        ggplot(taxon_mlst, aes(x="", y=x, fill=Label))+
        geom_bar(width=1, stat="identity", colour='white') +
        coord_polar("y", start=0) +
        scale_fill_manual(
            name="Scheme/ST",
            values=scale_values,
            labels=scale_labels
        ) +
        ggtitle(taxon) +
        theme_minimal() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold"),
            axis.text.x=element_blank()
        )
    print(p)
}
dev.off()
