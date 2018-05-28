#!/usr/bin/Rscript

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
# suppressMessages(library(superheat))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))

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

# test args
# args <- list(
#     stats='pg_stats.tsv',
#     meta_rds=c('sample_data.rds', 'sample_data_rest.rds'),
#     src_path="~/git_repos/Bacteria/SubProject_Pangenomes/scripts"
# )

## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--stats', help='', required=TRUE)
    parser$add_argument('--meta_rds', help='', required=TRUE, nargs='+')
    parser$add_argument('--obname', help='Output basename', required=TRUE)
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
data <- read_pg_stats(args$stats)
# Sample meta-data
tmp <- read_meta_rds(args$meta_rds)
meta <- tmp$meta
samples <- tmp$samples

## Plot data
# freq. matrix: taxon x freq.
freq <- data[,paste('freq', seq(10, 100, 10), sep='_')]
colnames(freq) <- sapply(colnames(freq), function(x){unlist(strsplit(x, '_'))[2]})

# re-shaped freq. matrix
freq2 <- melt.data.frame(cbind(Taxon=rownames(freq), freq), id.vars='Taxon')
colnames(freq2) <- c('Taxon', 'Freq', 'Num')
freq2$Freq <- as.numeric(as.character(freq2$Freq))
freq2$Text <- sapply(as.character(freq2$Taxon), function(x){ return(sprintf('%s (%d)', x, data[x,'SampleNum'])) })

# assembly/pan-genome size
tmp_median_asm  <- aggregate(meta[samples$main,'AssemblyLen'], by=list(Taxon=meta[samples$main,'TaxSpecies']), FUN=median); rownames(tmp_median_asm) <- tmp_median_asm$Taxon
tmp_mean_asm    <- aggregate(meta[samples$main,'AssemblyLen'], by=list(Taxon=meta[samples$main,'TaxSpecies']), FUN=mean); rownames(tmp_mean_asm) <- tmp_mean_asm$Taxon
tmp_min_asm     <- aggregate(meta[samples$main,'AssemblyLen'], by=list(Taxon=meta[samples$main,'TaxSpecies']), FUN=min); rownames(tmp_min_asm) <- tmp_min_asm$Taxon
tmp_max_asm     <- aggregate(meta[samples$main,'AssemblyLen'], by=list(Taxon=meta[samples$main,'TaxSpecies']), FUN=max); rownames(tmp_max_asm) <- tmp_max_asm$Taxon

asm_pg_size <- data.frame(
    Taxon=rownames(data),
    AssemblyMedian=tmp_median_asm[rownames(data),'x'],
    AssemblyMean=tmp_mean_asm[rownames(data),'x'],
    AssemblyMin=tmp_min_asm[rownames(data),'x'],
    AssemblyMax=tmp_max_asm[rownames(data),'x'],
    PangenomeSize=data[,'Total'],
    PangenomeCentroids=apply(freq, 1, sum)[rownames(data)],
    row.names=rownames(data), check.names=FALSE, stringsAsFactors=FALSE
)
# save
write.table(x=asm_pg_size, file=sprintf('%s.tsv', args$obname), row.names=FALSE, sep='\t', quote=FALSE)

## Plots
pt <-
    theme(
        axis.ticks=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=10, hjust=1, angle=90),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10, margin=margin(t=10, b=10)),
        axis.title.y=element_text(size=10, margin=margin(l=10, r=10)),
        strip.text.x=element_text(size=10),
        axis.text=element_text(colour="black"),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="#333333", fill=NA),
        strip.background=element_rect(colour="white", fill='white')
    )

pt2 <-
    theme(
        axis.ticks=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10, margin=margin(t=10, b=10)),
        axis.title.y=element_text(size=10, margin=margin(l=10, r=10)),
        strip.text.x=element_text(size=10),
        axis.text=element_text(colour="black"),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="#333333", fill=NA),
        strip.background=element_rect(colour="white", fill='white')
    )

#
Mb <- 1e+6
# Mb_trans <- trans_new(
#     name = 'Mb_trans',
#     transform = function(x){ x/Mb },
#     inverse = function(x){ x*Mb },
#     breaks = extended_breaks(),
#     format = format_format(),
#     domain = c(-Inf, Inf)
# )
asm_pg_size$AssemblyMedian  <- asm_pg_size$AssemblyMedian / Mb
asm_pg_size$AssemblyMean    <- asm_pg_size$AssemblyMean / Mb
asm_pg_size$AssemblyMin     <- asm_pg_size$AssemblyMin / Mb
asm_pg_size$AssemblyMax     <- asm_pg_size$AssemblyMax / Mb
asm_pg_size$PangenomeSize   <- asm_pg_size$PangenomeSize / Mb

p <-
    ggplot(data=asm_pg_size, aes(x=Taxon)) +
    # geom_label(aes(y=PangenomeSize, label=PangenomeCentroids), size=3, colour='black') +
    geom_point(aes(y=PangenomeSize), size=3, colour='black', shape=5) +
    geom_errorbar(aes(ymin=AssemblyMin, ymax=AssemblyMax), width=0.3, colour='black') +
    geom_point(aes(y=AssemblyMedian), size=2, shape=16) +
    scale_y_continuous(breaks=c(1,2,5,10,20,30), labels=c(1,2,5,10,20,30)) +
    # geom_point(aes(y=PangenomeSize, size=PangenomeCentroids), shape=23, fill='#CCCCCC', colour='#CCCCCC') +
    # scale_y_continuous( trans=Mb_trans ) +
    # scale_y_log10(
    #     breaks = trans_breaks("log10", function(x) 10^x),
    #     labels = trans_format("log10", math_format(10^.x))
    # ) +
    xlab('') + ylab('Assembly length (median, min/max) /\npan-genome length [Mb]') +
    # guides(size=guide_legend(title='Number of\ncentroids')) +
    pt2
    # pt

p2 <-
    ggplot(data=asm_pg_size, aes(x=Taxon, y=PangenomeCentroids)) +
    geom_bar(stat='identity') +
    geom_text(aes(y=PangenomeCentroids - 2500, label=PangenomeCentroids), size=3, colour='white') +
    # scale_y_log10(
    #     breaks = trans_breaks("log10", function(x) 10^x),
    #     labels = trans_format("log10", math_format(10^.x))
    # ) +
    xlab('') + ylab('# centroids') +
    # guides(size=guide_legend(title='Number of\ncentroids')) +
    pt

pdf(sprintf('%s.pdf', args$obname), height=9, width=9, onefile=FALSE)
# print(p)
grid.arrange(p, p2, ncol = 1)
dev.off()

# superheat plot
pdf(sprintf('%s_superheat.pdf', args$obname), height=6, width=9, onefile=FALSE)
p_sh <- superheat(
    X=log10(freq),
    # additional column plot
    # yt=log10(apply(freq, 2, sum)),
    # yt.axis.name="log10(# centroids)",
    # yt.plot.type="bar",
    # yt.bar.col="white",
    # yt.obs.col=rep("gray", ncol(freq)),
    # yt.plot.size=0.2,
    # additional row plot
    yr=data$SampleNum,
    yr.axis.name="# samples",
    yr.plot.type="bar",
    yr.bar.col="white",
    yr.obs.col=rep("gray", nrow(data)),
    yr.plot.size=0.2,
    # re-ordering and dendrograms
    pretty.order.rows=FALSE,
    pretty.order.cols=FALSE,
    row.dendrogram=FALSE,
    col.dendrogram=FALSE,
    # no smooth/scale
    smooth.heat=FALSE,
    scale=FALSE,
    # color palette
    heat.pal=c("#3399CC", '#FFFF99', '#FF6666'),
    heat.pal.values=c(0, 0.5, 1),
    heat.na.col='white',
    # grid
    grid.hline.size=0.3, grid.vline.size=0.3,
    grid.hline.col='#666666', grid.vline.col='#666666',
    # column labels
    # row labels
    left.label.text.size=4, left.label.size=0.5,
    # titles
    column.title='Centroid freq. [%]', row.title='Taxa', legend.title='log10(# centroids)',
)
dev.off()

# ggplot2
p <-
    ggplot(data=freq2, aes(x=Freq, y=Num)) +
    geom_bar(stat='identity') +
    facet_wrap(~Text, ncol=4) +
    scale_y_log10(
        breaks=trans_breaks("log10", function(x) 10^x),
        labels=trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
        breaks=seq(20, 100, 20)
    ) +
    ylab("Number of centroids") + xlab("Centroid frequency [%]") +
    pt
pdf(sprintf('%s_ggplot.pdf', args$obname), height=6, width=9, onefile=FALSE)
print(p)
dev.off()
