#!/usr/bin/Rscript

write("Rscript: Drug correlations\n", stdout())

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(reshape))
suppressMessages(library(Rtsne))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
# suppressMessages(library(superheat))
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
    parser$add_argument('--data', '-d', help='', required=TRUE)
    parser$add_argument('--obname', '-o', help='', required=TRUE)
    parser$add_argument('--min_num', help='', default=50, type="integer")
    parser$add_argument('--tsne_nopca', action="store_true")
    parser$add_argument('--tsne_perplex', default=10, type="integer")
    parser$add_argument('--tsne_theta', default=0, type="double")
    parser$add_argument('--dist', default='euclidean')
    parser$add_argument('--hclust', default='average')
    parser$add_argument('--src_path', help='')
    return(parser)
}

## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))
suppressMessages(source(sprintf('%s/utils/utils.R',args$src_path)))

## Log file
log_file <- sprintf('%s.log', args$obname)
write("Drug correlations\n", file=log_file, append=FALSE)
save_session_info(log_file)
save_args(args, log_file)

## Helper functions

prepare_data <- function(df, min_num=args$min_num, tax_rank='species'){
    # keep only isolates with taxon having at least N isolates
    write(sprintf('Preparing data: min. num. = %d, tax. rank %s', min_num, tax_rank), file=log_file, append=TRUE)
    # discard taxa
    df <- discard_taxa(df=df, min_num=min_num)
    # remove samples with taxon == "Rest"
    df <- df[df[,tax_rank] != "Rest",]
}

rm_const_cols <- function(df){
    # remove columns with only NAs or only one unique value
    return(
        df[,apply(df, 2, function(x){ length(unique(x)) > 1 }),drop=FALSE]
    )
}

write('Function to create binary resistance profiles:', file=log_file, append=TRUE)
res2bin <- function(df, one_values=c('R')){
    return(apply(df, c(1,2), function(x){ ifelse(is.na(x), NA, ifelse(x %in% one_values, 1, 0)) }))
}
write(capture.output(print(res2bin)), file=log_file, append=TRUE)

superheat_cor <- function(cor_mat, dist_method=args$dist, title='', draw_legend=FALSE, draw_text=FALSE){
    # complete correlation range - color transitions
    complete_pal_values <- seq(-1, 1, 0.1)
    # palette for above values
    complete_pal <- colorRampPalette(rev(c("#3399CC", 'white', '#FF6666')))(length(complete_pal_values))
    # round correlation values, i.e. can match to above palette values
    cor_mat <- round(cor_mat, 1)
    cor_mat_pal <- (complete_pal_values >= min(cor_mat) & complete_pal_values <= max(cor_mat))
    # drug class colors
    drug_colors <- drug_classes[tabs$drugs[colnames(cor_mat),'Class'],'Color']
    # plot
    superheat(
        X=cor_mat,
        X.text=apply(cor_mat, c(1,2), function(x){ ifelse(draw_text, sprintf('%.1f', x), '') }),
        X.text.size=2,
        # re-order rows/columns and dendrograms
        pretty.order.rows=TRUE, pretty.order.cols=TRUE,
        row.dendrogram=FALSE,
        col.dendrogram=FALSE,
        clustering.method='hierarchical', dist.method=dist_method, linkage.method="average",
        #
        smooth.heat=FALSE,
        scale=FALSE,
        # color palette
        legend=draw_legend,
        heat.pal=complete_pal[cor_mat_pal],
        heat.pal.values=scales::rescale(x=complete_pal_values[cor_mat_pal], to=c(0, 1)),
        # grid
        grid.hline.size=0.3, grid.vline.size=0.3,
        grid.hline.col='#666666', grid.vline.col='#666666',
        # column labels
        bottom.label="variable", bottom.label.col=drug_colors, bottom.label.text.size=2, bottom.label.size=0.075,
        # row labels
        left.label="variable", left.label.col=drug_colors, left.label.text.size=2, left.label.size=0.075,
        title=title, title.size=11,
        padding=0.1
    )
}

##
log_hrule(log_file)

## Data
tabs <- split_data(tab=read_data(file_name=args$data))

## Samples
# any location and with MICs
ids_all <-  rownames(tabs$mics)[sapply(rownames(tabs$mics), function(s_id){
    !all(is.na(tabs$mics[s_id,])) # has MIC profile
})]

# log
write(sprintf("Filtering by MIC profile: %d remained", length(ids_all)), file=log_file, append=TRUE)

## Data for correlation computation
df_mics <- data.frame(
    species=convert_taxa(tabs$tax[ids_all,'species']), # converted taxa names
    tabs$mics[ids_all,], # MIC values
    check.names=FALSE, stringsAsFactors=FALSE, row.names=ids_all
)
df_mics <- prepare_data(df=df_mics, min_num=args$min_num, tax_rank='species') # filter isolates by taxon
# log
write(sprintf(
    "Filtering by number of samples per taxon: %d samples, %d taxa",
    nrow(df_mics), length(unique(df_mics$species))
    ),
    file=log_file, append=TRUE
)
ids_all <- rownames(df_mics) # because some were removed (taxonomy == "Rest")

# checks
testit::assert(all(sapply(ids_all, function(x){x %in% rownames(df_mics)})))
write(sprintf("Final: %d samples", length(ids_all)), file=log_file, append=TRUE)

# all unique taxa
taxa <- sort(unique(df_mics$species))

## Correlation matrices
cor_mics <- list(
    All=cor(x=rm_const_cols(df_mics[ids_all, -1]), method='spearman')
)

for(taxon in taxa){
    # get IDs
    taxon_ids_all <- ids_all[ df_mics[ids_all,'species'] == taxon ]
    # check
    testit::assert(unique(df_mics[taxon_ids_all, 'species'])==taxon)
    # create value tables
    taxon_mics_all <- rm_const_cols(df_mics[taxon_ids_all, -1])
    # compute correlation
    taxon_cor_mics_all <- cor(x=taxon_mics_all, method='spearman')
    # append to lists
    cor_mics <- c(cor_mics, list(list(All=taxon_cor_mics_all)))
    names(cor_mics)[length(cor_mics)] <- taxon
}

## Stat.s
## AZT - CPE/CAX/CFT/CAZ vs. rest
log_hrule(log_file)
write("Special case: AZT: median cor. CPE/CAX/CFT/CAZ vs. median cor. rest\\AZT", file=log_file, append=TRUE)
azt_gr1 <- c('CPE','CAX','CFT','CAZ')
azt_gr2 <- setdiff(tabs$drug$Abbr, c('AZT', azt_gr1))

write(sprintf(
    "\t%.1f vs. %.1f: All taxa",
    median(cor_mics$All['AZT',azt_gr1]),
    median(cor_mics$All['AZT',intersect(colnames(cor_mics$All),azt_gr2)])
), file=log_file, append=TRUE)
for(taxon in taxa){
    if(! 'AZT' %in% colnames(cor_mics[[taxon]]$All)){
        write(sprintf('\t NA vs. NA: %s (AZT was removed)', taxon), file=log_file, append=TRUE)
        next
    }
    write(sprintf(
        "\t%.1f vs. %.1f: %s",
        median(cor_mics[[taxon]]$All['AZT',azt_gr1]),
        median(cor_mics[[taxon]]$All['AZT',intersect(colnames(cor_mics[[taxon]]$All),azt_gr2)]),
        taxon
    ), file=log_file, append=TRUE)
}

## Heatmaps
log_hrule(log_file)
write("MIC profile based", file=log_file, append=TRUE)

# color key
pdf(sprintf('%s_mics_superheat_colorkey.pdf', args$obname), width=16, height=1.5)
complete_pal_values <- seq(-1, 1, 0.1)
complete_pal <- colorRampPalette(rev(c("#3399CC", 'white', '#FF6666')))(length(complete_pal_values))
df <- data.frame(y=rep(1, length(complete_pal_values)),x=complete_pal_values)
p_kc <-
    ggplot(data=df, aes(y=y, x=x, fill=x)) +
    geom_tile() +
    scale_fill_gradientn(colours=complete_pal, guide=FALSE) +
    xlab('') + ylab('') +
    theme_bw() +
    theme(
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(size=28)
    )
print(p_kc)
dev.off()

# superheat
pdf(sprintf('%s_mics_superheat.pdf', args$obname), width=6, height=7)
# All
# superheat_cor(cor_mics$All, args$dist, 'Drug correlation\nAll taxa')
superheat_cor(cor_mics$All, args$dist, 'All taxa')
# Per species
for(taxon in taxa){
    if(any(dim(cor_mics[[taxon]]$All) < 10)){
        write(sprintf("MICs: Skipped %s: %s", taxon, paste(dim(cor_mics[[taxon]]$All), collapse=' x ')), file=log_file, append=TRUE)
        next
    }
    # superheat_cor(cor_mics[[taxon]]$All, args$dist, sprintf('Drug correlation\n%s', taxon))
    superheat_cor(cor_mics[[taxon]]$All, args$dist, gsub(' ' ,'\n', sprintf('%s', taxon)))
}
dev.off()

pdf(sprintf('%s_mics_superheat_with_legend.pdf', args$obname), width=6, height=7)
# All
superheat_cor(cor_mics$All, args$dist, 'All taxa', draw_legend=TRUE, draw_text=TRUE)
# Per species
for(taxon in taxa){
    if(any(dim(cor_mics[[taxon]]$All) < 10)){
        write(sprintf("MICs: Skipped %s: %s", taxon, paste(dim(cor_mics[[taxon]]$All), collapse=' x ')), file=log_file, append=TRUE)
        next
    }
    superheat_cor(cor_mics[[taxon]]$All, args$dist, gsub(' ' ,'\n', sprintf('%s', taxon)), draw_legend=TRUE, draw_text=TRUE)
}
dev.off()

# tiff_count <- 1
# tiff(filename=sprintf('%s/tiff/%s_mics_superheat_%d.tiff', dirname(args$obname), basename(args$obname), tiff_count), width=6, height=7, res=300, units = 'in', compression="lzw")
# superheat_cor(cor_mics$All, args$dist, 'Drug correlation\nAll taxa')
# # superheat_cor(cor_mics$All, args$dist, '')
# dev.off()
# tiff_count <- tiff_count + 1
#
# for(taxon in taxa){
#     if(any(dim(cor_mics[[taxon]]$All) < 10)){
#         next
#     }
#     tiff(filename=sprintf('%s/tiff/%s_mics_superheat_%d.tiff', dirname(args$obname), basename(args$obname), tiff_count), width=6, height=7, res=300, units = 'in', compression="lzw")
#     # superheat_cor(cor_mics[[taxon]]$All, args$dist, sprintf('Drug correlation\n%s', taxon))
#     superheat_cor(cor_mics[[taxon]]$All, args$dist, sprintf('%s', taxon))
#     dev.off()
#     tiff_count <- tiff_count + 1
# }
