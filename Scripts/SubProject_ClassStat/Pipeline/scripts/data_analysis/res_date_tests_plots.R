#!/usr/bin/Rscript

write("Rscript: resistance-date test plots\n", stdout())

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(reshape2))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
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
    parser$add_argument('--rds', required=TRUE)
    parser$add_argument('--obname', required=TRUE)
    parser$add_argument('--src_path')
    parser$add_argument('--alpha', default=0.05, type="double")
    return(parser)
}

## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## Help functions
suppressMessages(source(sprintf('%s/utils/utils.R',args$src_path)))

## Log file
log_file <- sprintf('%s.log', args$obname)
write("resistance-date tests\n", file=log_file, append=FALSE)
save_session_info(log_file)
save_args(args, log_file)

## Data
test_results <- readRDS(args$rds)

## P-values
pvs_adj_log <- -log10(test_results$pvalues_adj)
pvs_adj_log_vec <- unlist(pvs_adj_log)
pvs_adj_log[pvs_adj_log == Inf] <- max(pvs_adj_log_vec[is.finite(pvs_adj_log_vec)], na.rm=TRUE)
pvs_adj_log[is.na(pvs_adj_log)] <- 0 # won't affect dist. for clustering
# heat_pal_values <- c(
#     0, # no test
#     -log10(args$alpha), # > is sign.
#     quantile(pvs_adj_log[pvs_adj_log > -log10(args$alpha)], c(0, 0.5, 0.75, 1))
# )

## Heatmap text
hm_text <- matrix('', nrow=nrow(test_results$pvalues_adj), ncol=ncol(test_results$pvalues_adj), dimnames=dimnames(test_results$pvalues_adj))
for(taxon in rownames(hm_text)){
    for(drug in colnames(hm_text)){
        if(is.na(test_results$pvalues_adj[taxon, drug])){
            next
        }
        if(test_results$pvalues_adj[taxon, drug] < args$alpha){
            if( test_results$tests[[sprintf('%s:%s', taxon, drug)]]$estimate < 0 ){
                hm_text[taxon, drug] <- 'x' #sprintf('\u2191')
            } else {
                hm_text[taxon, drug] <- 'o' #sprintf('\u2193')
            }
        }
    }
}

## Heatmap: taxa x drugs
pdf(sprintf("%s_heatmap.pdf", args$obname), width=11, height=9, onefile=FALSE)
superheat(
    # plot data
    X=pvs_adj_log,
    # X.text=apply(test_results$pvalues_adj, c(1,2), function(x){ ifelse(is.na(x), '', sprintf('%.0e', x)) }),
    X.text=hm_text,
    # yt=NULL,
    # yr=NULL,
    # dendrograms/re-ordering
    # membership.rows=NULL,
    # membership.cols=NULL,
    pretty.order.rows=TRUE, pretty.order.cols=TRUE,
    row.dendrogram=TRUE, col.dendrogram=TRUE,
    # n.clusters.rows=NULL, n.clusters.cols=NULL,
    clustering.method="hierarchical", #c("kmeans", "hierarchical"),
    dist.method="euclidean", # c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
    linkage.method="average",
    # order.cols=NULL, order.rows=NULL,
    # scaling
    smooth.heat=FALSE,
    scale=FALSE,
    # axis labels
    left.label="variable", bottom.label="variable",
    # color scheme
    heat.col.scheme="viridis", # c("viridis", "red", "purple", "blue", "grey", "green"
    # heat.pal=c("white", "gray50", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"),
    # heat.pal.values=rescale(x=heat_pal_values, to=c(0, 1)),
    heat.na.col="white",
    heat.lim=c(-log10(args$alpha), max(pvs_adj_log, na.rm=TRUE)),
    # text in the cells
    # X.text.size=2, X.text.col="black", X.text.angle=0,
    # legend
    legend=TRUE, legend.height=0.1, legend.width=1.5, legend.text.size=12,
    # grid
    grid.hline=TRUE, grid.vline=TRUE,
    # grid.hline.size=0.5, grid.vline.size=0.5,
    # grid.hline.col="black", grid.vline.col="black",
    # force.grid.hline=F, force.grid.vline=F,
    # smooth
    # smoothing.method=c("loess", "lm"), smooth.se=TRUE,
    # sub-plot axes
    # yt.plot.type=c("scatter", "bar", "boxplot", "scattersmooth", "smooth", "scatterline", "line"),
    # yr.plot.type=c("scatter", "bar", "boxplot", "scattersmooth", "smooth", "scatterline", "line"),
    # yt.axis=T, yr.axis=T,
    # yt.num.ticks=3, yr.num.ticks=3,
    # yt.plot.size=0.3, yr.plot.size=0.3,
    # yt.axis.name=NULL, yr.axis.name=NULL,
    # yr.axis.size=10, yt.axis.size=10,
    # yr.axis.name.size=10, yt.axis.name.size=10,
    # yr.axis.name.angle=NULL, yt.axis.name.angle=NULL,
    # yt.obs.col=NULL, yr.obs.col=NULL,
    # yt.cluster.col=NULL, yr.cluster.col=NULL,
    # yt.bar.col=NULL, yr.bar.col=NULL,
    # yt.point.size=2, yt.point.alpha=1, yr.point.size=2, yr.point.alpha=1,
    # yr.line.col=NULL, yt.line.col=NULL, yr.line.size=NULL, yt.line.size=NULL,
    # axes
    # bottom.label.text.size=5, left.label.text.size=5,
    bottom.label.text.angle=90, left.label.text.angle=NULL,
    bottom.label.size=0.1, left.label.size=0.4,
    # left.label.col=NULL, bottom.label.col=NULL,
    # left.label.text.col=NULL, bottom.label.text.col=NULL,
    left.label.text.alignment="center", bottom.label.text.alignment="center", # c("center", "left", "right")
    # force.left.label=F, force.bottom.label=F,
    # titles
    padding=1,
    column.title="-log10(adj. p-value)", row.title="",
    column.title.size=5, row.title.size=5,
    title=NULL, title.size=5,
    #
    print.plot=TRUE
)
dev.off()

pdf(sprintf("%s_heatmap_text.pdf", args$obname), width=11, height=9, onefile=FALSE)
superheat(
    # plot data
    X=pvs_adj_log,
    X.text=apply(test_results$pvalues_adj, c(1,2), function(x){ ifelse(is.na(x), '', sprintf('%.0e', x)) }),
    # yt=NULL,
    # yr=NULL,
    # dendrograms/re-ordering
    # membership.rows=NULL,
    # membership.cols=NULL,
    pretty.order.rows=TRUE, pretty.order.cols=TRUE,
    row.dendrogram=TRUE, col.dendrogram=TRUE,
    # n.clusters.rows=NULL, n.clusters.cols=NULL,
    clustering.method="hierarchical", #c("kmeans", "hierarchical"),
    dist.method="euclidean", # c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
    linkage.method="average",
    # order.cols=NULL, order.rows=NULL,
    # scaling
    smooth.heat=FALSE,
    scale=FALSE,
    # axis labels
    left.label="variable", bottom.label="variable",
    # color scheme
    heat.col.scheme="viridis", # c("viridis", "red", "purple", "blue", "grey", "green"
    # heat.pal=c("white", "gray50", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"),
    # heat.pal.values=rescale(x=heat_pal_values, to=c(0, 1)),
    heat.na.col="white",
    heat.lim=c(-log10(args$alpha), max(pvs_adj_log, na.rm=TRUE)),
    # text in the cells
    X.text.size=2, X.text.col="black", X.text.angle=0,
    # legend
    legend=TRUE, legend.height=0.1, legend.width=1.5, legend.text.size=12,
    # grid
    grid.hline=TRUE, grid.vline=TRUE,
    # grid.hline.size=0.5, grid.vline.size=0.5,
    # grid.hline.col="black", grid.vline.col="black",
    # force.grid.hline=F, force.grid.vline=F,
    # smooth
    # smoothing.method=c("loess", "lm"), smooth.se=TRUE,
    # sub-plot axes
    # yt.plot.type=c("scatter", "bar", "boxplot", "scattersmooth", "smooth", "scatterline", "line"),
    # yr.plot.type=c("scatter", "bar", "boxplot", "scattersmooth", "smooth", "scatterline", "line"),
    # yt.axis=T, yr.axis=T,
    # yt.num.ticks=3, yr.num.ticks=3,
    # yt.plot.size=0.3, yr.plot.size=0.3,
    # yt.axis.name=NULL, yr.axis.name=NULL,
    # yr.axis.size=10, yt.axis.size=10,
    # yr.axis.name.size=10, yt.axis.name.size=10,
    # yr.axis.name.angle=NULL, yt.axis.name.angle=NULL,
    # yt.obs.col=NULL, yr.obs.col=NULL,
    # yt.cluster.col=NULL, yr.cluster.col=NULL,
    # yt.bar.col=NULL, yr.bar.col=NULL,
    # yt.point.size=2, yt.point.alpha=1, yr.point.size=2, yr.point.alpha=1,
    # yr.line.col=NULL, yt.line.col=NULL, yr.line.size=NULL, yt.line.size=NULL,
    # axes
    # bottom.label.text.size=5, left.label.text.size=5,
    bottom.label.text.angle=90, left.label.text.angle=NULL,
    bottom.label.size=0.1, left.label.size=0.4,
    # left.label.col=NULL, bottom.label.col=NULL,
    # left.label.text.col=NULL, bottom.label.text.col=NULL,
    left.label.text.alignment="center", bottom.label.text.alignment="center", # c("center", "left", "right")
    # force.left.label=F, force.bottom.label=F,
    # titles
    padding=1,
    column.title="-log10(adj. p-value)", row.title="",
    column.title.size=5, row.title.size=5,
    title=NULL, title.size=5,
    #
    print.plot=TRUE
)
dev.off()

## Boxplots
box_data <- test_results$test_data # get data
box_data <- box_data[, setdiff(colnames(box_data), c('Date'))] # rm date column
box_data <- reshape2::melt(box_data, id.vars=c('Year', 'species')) # re-shape: all drug in one column
for(taxon in rownames(test_results$pvalues_adj)){ # filter by taxon and drug: keep only sign.
    for(drug in colnames(test_results$pvalues_adj)){
        rm_taxon_drug <- FALSE
        if(is.na(test_results$pvalues_adj[taxon, drug])){
            rm_taxon_drug <- TRUE
        } else if (test_results$pvalues_adj[taxon, drug] >= args$alpha){
            rm_taxon_drug <- TRUE
        }
        if(rm_taxon_drug){
            box_data <- box_data[!(box_data$species==taxon & box_data$variable==drug),]
        }
    }
}
box_data$value <- sapply(box_data$value, function(x){ ifelse(x=='R', 'R', 'S+I') }) # modify values
p <-ggplot(data=box_data, aes(x=value, y=Year, fill=value)) +
    scale_fill_manual(values=c("R"="#FF6666", "S+I"="#3399FF")) +
    geom_boxplot(colour='black') +
    facet_grid(species~variable) +
    xlab('Resistance value') + ylab('Isolate collection year') +
    # ggtitle(taxon) +
    theme_bw()
pdf(sprintf("%s_boxplots.pdf", args$obname), width=22, height=20)
print(p)
dev.off()

p <-ggplot(data=box_data, aes(x=Year, colour=value)) +
    scale_colour_manual(values=c("R"="#FF6666", "S+I"="#3399FF")) +
    geom_density(fill=NA) +
    facet_grid(species~variable) +
    xlab('Resistance value') + ylab('Isolate collection year') +
    # ggtitle(taxon) +
    theme_bw()
pdf(sprintf("%s_density.pdf", args$obname), width=22, height=20)
print(p)
dev.off()
