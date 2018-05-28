#!/usr/bin/Rscript

min_prs_pct <- 90

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
# suppressMessages(library(superheat))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))

suppressMessages(library(magrittr))
source('~/R/custom/superheat-0.1.0/R/addplot.R')
source('~/R/custom/superheat-0.1.0/R/cluster.R')
source('~/R/custom/superheat-0.1.0/R/heatmap.R')
source('~/R/custom/superheat-0.1.0/R/labels.R')
source('~/R/custom/superheat-0.1.0/R/layout.R')
source('~/R/custom/superheat-0.1.0/R/plot-types.R')
source('~/R/custom/superheat-0.1.0/R/superheat.R')
source('~/R/custom/superheat-0.1.0/R/themes.R')
source('~/R/custom/superheat-0.1.0/R/utils.R')

## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--hits_rds', help='', required=TRUE)
    parser$add_argument('--annot', help='', required=TRUE)
    parser$add_argument('--meta_rds', help='', required=TRUE, nargs='+')
    parser$add_argument('--obname', help='', required=TRUE)
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
objs <- readRDS(args$hits_rds)
# resfams annotations
annot <- read.csv(file=args$annot, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
rownames(annot) <- annot[,'Resfam ID']
# sample meta-data
tmp <- read_meta_rds(args$meta_rds)
meta <- tmp$meta
samples <- tmp$samples
write(sprintf('Meta data: %d samples (%d, rest %d)', nrow(meta), length(samples$main), length(samples$rest)), log_file, append=TRUE)

write(capture.output(print(table(round(objs$ov_g_s_prs)))), log_file, append=TRUE)

## NOTE Stat.s
# for all species: how many with mean count >= 1
for(taxon in rownames(objs$ov_g_sp)){
    testit::assert(ncol(objs$ov_g_sp) == nrow(annot))
    taxon_num <- sum(objs$ov_g_sp[taxon,] >= 1)
    taxon_pct <- 100 * taxon_num / ncol(objs$ov_g_sp)
    write(sprintf('\n%s: %.1f (%d of %d) of Resfams genes', taxon, taxon_pct, taxon_num, ncol(objs$ov_g_sp)), log_file, append=TRUE)
}

# S. aureus: how many with mean count >= 1
# testit::assert(ncol(objs$ov_g_sp) == nrow(annot))
# saureus_num <- sum(objs$ov_g_sp['Staphylococcus aureus',] >= 1)
# saureus_pct <- 100 * saureus_num / ncol(objs$ov_g_sp)
# write(sprintf('\nS. aureus: %.1f (%d of %d) of Resfams genes', saureus_pct, saureus_num, ncol(objs$ov_g_sp)), log_file, append=TRUE)
# # S. aureus: mean count > 60
# write(sprintf('\nS. aureus: mean hit count > 60: %s', paste(colnames(objs$ov_g_sp)[objs$ov_g_sp['Staphylococcus aureus',] > 60], collapse=',')), log_file, append=TRUE)
# Mean cout per species in RF0115
write(sprintf(
    '\nMean hit count for RF0115 per species:\n\t%s',
    paste(sapply(rownames(objs$ov_g_sp), function(x){ sprintf('%s: %.1f', x, objs$ov_g_sp[x,'RF0115']) }), collapse='\n\t')
    ), log_file, append=TRUE
)
# In S. aureus only
in_saureus_only <-
    apply(objs$ov_g_sp[setdiff(rownames(objs$ov_g_sp), 'Staphylococcus aureus'),], 2, function(x){ mean(x) < 1 }) &
    objs$ov_g_sp['Staphylococcus aureus',] >= 1
write(sprintf(
    '\nIn S. aureus only (mean hit count >= 1, mean(other species) < 1): %s',
    paste(colnames(objs$ov_g_sp)[in_saureus_only], collapse=',')
    ), log_file, append=TRUE
)
# Resistance genes found in at least 90% of isolates
in_ge90 <- names(objs$ov_g_s_prs)[objs$ov_g_s_prs >= min_prs_pct]
write(sprintf(
    '\nResfams genes found in >= %dpct of isolates: %d: %s',
    min_prs_pct, length(in_ge90), paste(in_ge90, collapse=',')
    ), log_file, append=TRUE
)

## NOTE Plots
write('\nCreate plots', log_file, append=TRUE)
pdf(sprintf('%s.pdf', args$obname), width=9, height=6)

## Plot: Samples x Resfans
plot_ov_g_s_ss <- function(min_prs=0, max_prs=101){
    # which columns to plot (names)
    ov_g_s_ss <- names(objs$ov_g_s_prs)[objs$ov_g_s_prs >= min_prs & objs$ov_g_s_prs < max_prs]
    # taxa for rows
    ov_g_s_taxa <- abbr_taxa(convert_taxa(meta[rownames(objs$ov_g_s),'TaxSpecies']))
    # unique taxa
    u_taxa <- sort(unique(ov_g_s_taxa))
    # empty string (for plotting taxa)
    s <- paste(rep(' ', 25), collapse='')
    # modify taxa for plotting (so they do not overlap)
    ov_g_s_taxa <- sapply(ov_g_s_taxa, function(x){ ifelse((which(u_taxa==x) %% 2) == 0, paste(x, s, sep=''), paste(s, x, sep='')) })
    # custom values for color breaks
    my_pal_values <- scales::rescale(c(0,0.5,0.75,1,max(objs$ov_g_s[,ov_g_s_ss])), to=c(0,1))

    superheat(
        X=objs$ov_g_s[,ov_g_s_ss],
        membership.rows=ov_g_s_taxa,
        pretty.order.rows=FALSE,
        pretty.order.cols=TRUE,
        row.dendrogram=FALSE,
        col.dendrogram=FALSE,
        smooth.heat=FALSE,
        scale=FALSE,
        linkage.method="average",
        yt = round(objs$ov_g_s_prs[ov_g_s_ss]),
        yt.axis.name = "Presence [%]",
        yt.plot.type = "bar",
        heat.col.scheme='viridis',
        heat.pal.values=my_pal_values,
        legend.breaks=c(0, 1, 2, 5, 10),
        legend=TRUE,
        grid.hline.size=0.1, grid.vline.size=0.1,
        bottom.label='variable', bottom.label.text.angle='90', bottom.label.text.size=2, bottom.label.size=0.2,
        left.label.text.size=3, left.label.size=0.3,
        row.title='Samples', column.title='Resfams core'
    )
}
plot_ov_g_s_ss(0,10)
plot_ov_g_s_ss(10,15)
plot_ov_g_s_ss(15,50)
plot_ov_g_s_ss(50,75)
plot_ov_g_s_ss(75,101)

## Plot: Species x Resfams
my_pal_values <- scales::rescale(c(0,0.5,0.75,1,max(objs$ov_g_sp)), to=c(0,1))
superheat(
    X=objs$ov_g_sp,
    pretty.order.rows=TRUE,
    pretty.order.cols=TRUE,
    row.dendrogram=TRUE,
    col.dendrogram=FALSE,
    smooth.heat=FALSE,
    scale=FALSE,
    linkage.method="average",
    legend=TRUE,
    yt = round(objs$ov_g_s_prs[colnames(objs$ov_g_sp)]),
    yt.obs.col = sapply(objs$ov_g_s_prs[colnames(objs$ov_g_sp)], function(x){ ifelse(x >= min_prs_pct, '#FF3333', '#333333') }),
    yt.axis.name = "Presence [%]",
    yt.plot.type = "bar",
    heat.col.scheme='viridis',
    heat.pal.values=my_pal_values,
    legend.breaks=c(0, 1, 2, 5, 10),
    grid.hline.size=0.1, grid.vline.size=0.1,
    left.label.text.size=3, left.label.size=0.5,
    row.title='Taxa', column.title='Resfams core'
)

dev.off()
