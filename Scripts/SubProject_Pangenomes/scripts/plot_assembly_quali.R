#!/usr/bin/Rscript

CUTOFF_CONTIG <- 1000
CUTOFF_L50    <- 200
CUTOFF_N50    <- 5000

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))

# args <- list(
#     meta_rds=c('sample_data.rds', 'sample_data_rest.rds'),
#     src_path="~/git_repos/Bacteria/SubProject_Pangenomes/scripts"
# )


## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--meta_rds', help='', required=TRUE, nargs='+')
    parser$add_argument('--asm_cov', help='', required=TRUE)
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
tmp <- read_meta_rds(args$meta_rds)
meta <- tmp$meta
samples <- tmp$samples
asm_cov <- read.csv(file=args$asm_cov, header=FALSE, stringsAsFactors=FALSE, sep='\t'); colnames(asm_cov) <- c('sample', 'mean_cov'); rownames(asm_cov) <- asm_cov$sample

tmp_passed <- meta$AssemblyNum <= CUTOFF_CONTIG & meta$AssemblyL50 <= CUTOFF_L50 & meta$AssemblyN50 >= CUTOFF_N50
write(sprintf('Passing quality criteria: %d (%.1f pct of %d, %.1f pct of %d w/ assembly)',
    sum(tmp_passed, na.rm=TRUE),
    100 * sum(tmp_passed, na.rm=TRUE) / nrow(meta), nrow(meta),
    100 * sum(tmp_passed, na.rm=TRUE) / sum(!is.na(meta$AssemblyNum)), sum(!is.na(meta$AssemblyNum))
), log_file, append=TRUE)

## Main taxa
# from sampels used in pan-genomes get unique species taxa
main_taxa <- convert_taxa(sort(unique(meta[samples$main, 'TaxSpecies'])))
write(sprintf('Main species taxa: %s', paste(main_taxa, collapse=', ')), log_file, append=TRUE)
# Lab strain samples from S. aureus dataset
lab_strain_samples <- rownames(meta)[meta$DataSet=='saureus' & meta$Source=='Lab strain']
write(sprintf('S. aureus lab strain samples: %s', paste(lab_strain_samples, collapse=', ')), log_file, append=TRUE)

## Plot data
df <- data.frame(
    # sample ID
    ID=rownames(meta),
    # Assembly stat.s
    meta[,grepl('Assembly', colnames(meta))],
    # Mean assembly coverage
    AssemblyMeanCov=asm_cov[rownames(meta),'mean_cov'],
    # Taxon
    Taxon=convert_taxa(meta$TaxSpecies),
    #
    row.names=rownames(meta), stringsAsFactors=FALSE, check.names=FALSE
)
# rm samples without assembly
df <- df[!is.na(df$AssemblyNum),]
write(sprintf(
    'Removed samples w/o assembly, kept %d (%.1f pct of %d)',
    nrow(df), 100 * nrow(df) / nrow(meta), nrow(meta)
    ), log_file, append=TRUE
)
testit::assert(length(intersect(df$ID, asm_cov$sample)) == length(df$ID))
# discard taxa
df$Taxon <- sapply(df$Taxon, function(x){ ifelse(x %in% main_taxa, x, 'Other') })
# abbr. taxon
# df$Taxon <- abbr_taxa(df$Taxon)
# Samples used in pan-genomes
df$InPG <- sapply(df$ID, function(x){ ifelse(x %in% samples$main, 'Yes', 'No') })

# Re-shape plot data
df_melt <- melt.data.frame(
    data=df[,setdiff(colnames(df), c('AssemblyLen', 'AssemblyMaxLen', 'AssemblyNs'))],
    id.vars=c('ID', 'Taxon', 'InPG')
)
# variable -> change names
# df_melt$variable <- sapply(as.character(df_melt$variable), function(x){
#     switch(x,
#         'AssemblyGC' ={ return('GC content [%]') },
#         'AssemblyNum'={ return('Number of contigs') },
#         'AssemblyN50'={ return('N50 [bp]') },
#         'AssemblyL50'={ return('L50') },
#         'AssemblyMeanCov'={ return('Mean assembly coverage') },
#         { return(x) }
#     )
# })
# order variable and taxon
df_melt$variable <- factor(
    x=as.character(df_melt$variable),
    levels=c('AssemblyGC','AssemblyMeanCov','AssemblyNum','AssemblyL50','AssemblyN50'),
    ordered=TRUE
)
df_melt$Taxon <- factor(
    x=as.character(df_melt$Taxon),
    levels=c(sort(setdiff(unique(df_melt$Taxon), 'Other')), 'Other'),
    ordered=TRUE
)

# threshold values
df_th <- data.frame(
    variable=factor(
        x=unique(as.character(df_melt$variable)),
        levels=c('AssemblyGC','AssemblyMeanCov','AssemblyNum','AssemblyL50','AssemblyN50'),
        ordered=TRUE
    ),
    Y=NA
)
rownames(df_th) <- as.character(df_th$variable)
df_th['AssemblyNum', 'Y'] <- CUTOFF_CONTIG
df_th['AssemblyL50', 'Y'] <- CUTOFF_L50
df_th['AssemblyN50', 'Y'] <- CUTOFF_N50

# threshold rect
df_rect <- data.frame(
    variable=factor(
        x=unique(as.character(df_melt$variable)),
        levels=c('AssemblyGC','AssemblyMeanCov','AssemblyNum','AssemblyL50','AssemblyN50'),
        ordered=TRUE
    ),
    stringsAsFactors=FALSE
)
rownames(df_rect) <- df_rect$variable
df_rect['AssemblyNum', 'ymax'] <- df_th['AssemblyNum','Y']
df_rect['AssemblyL50', 'ymax'] <- df_th['AssemblyL50','Y']
df_rect['AssemblyN50', 'ymax'] <- Inf
df_rect['AssemblyNum', 'ymin'] <- 0
df_rect['AssemblyL50', 'ymin'] <- 0
df_rect['AssemblyN50', 'ymin'] <- df_th['AssemblyN50','Y']

pt <-
    theme(
        axis.ticks=element_blank(),
        panel.grid.major=element_line(colour="#CCCCCC"),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=9, hjust=1, angle=90),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10, margin=margin(t=10, b=10), face='bold'),
        axis.title.y=element_text(size=10, margin=margin(l=10, r=10)),
        strip.text.x=element_text(size=10),
        axis.text=element_text(colour="black"),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="#333333", fill=NA),
        strip.background=element_rect(colour="white", fill='white')
    )

s <- 'isolates passing threshold'
passed_cn <- sum(df$AssemblyNum <= CUTOFF_CONTIG)
passed_lf <- sum(df$AssemblyL50 <= CUTOFF_L50)
passed_nf <- sum(df$AssemblyN50 >= CUTOFF_N50)
passed_all<- sum(
    df$AssemblyNum <= CUTOFF_CONTIG &
    df$AssemblyL50 <= CUTOFF_L50 &
    df$AssemblyN50 >= CUTOFF_N50
)
total_num <- length(unique(df_melt$ID))

custom_labeller <- c(
    'AssemblyGC'     = 'GC content [%]',
    'AssemblyNum'    = sprintf('Number of contigs\n(cut-off value = %d)', CUTOFF_CONTIG),
    'AssemblyN50'    = sprintf('N50 [bp]\n(cut-off value = %d bp)', CUTOFF_N50),
    'AssemblyL50'    = sprintf('L50\n(cut-off value = %d)', CUTOFF_L50),
    'AssemblyMeanCov'= 'Mean assembly coverage'
)

p1 <- # all samples
    ggplot(data=df_melt, aes(x=Taxon, y=value)) + # ggplot(data=df_melt[df_melt$Taxon != 'Other',], aes(x=Taxon, y=value)) +
    geom_rect(data=df_rect, inherit.aes=FALSE, aes(ymin=ymin, ymax=ymax), xmin=-Inf, xmax=Inf, alpha=0.5, fill='#33CC00') +
    geom_boxplot() +
    geom_hline(data=df_th, aes(yintercept=Y), linetype="dotted") +
    facet_wrap(~variable, ncol=1, scales='free_y', labeller=labeller(variable=custom_labeller)) +
    scale_y_sqrt(label=comma) +
    ylab('') + xlab('')

p2 <- # only samples in pan-genomes
    ggplot(data=df_melt[df_melt$InPG == 'Yes',], aes(x=Taxon, y=value)) +
    geom_boxplot() +
    geom_hline(data=df_th, aes(yintercept=Y), linetype="dotted") +
    facet_wrap(~variable, ncol=1, scales='free_y', labeller=labeller(variable=custom_labeller)) +
    xlab('') + ylab('')

pdf(sprintf('%s.pdf', args$obname), width=9, height=9)
ggplot_theme_setting <- theme_set(theme_bw() + pt)
plot_grid(p1, NULL, align='h', ncol=2, nrow=1, rel_widths=c(0.7, 0.3)) +
    draw_text(
        sprintf('%.1f%%', 100 * passed_cn / total_num),
        x = 0.84, y = 0.57, size = 28,
        colour='#33CC00', alpha=0.75
    ) +
    draw_text(
        sprintf('%.1f%%', 100 * passed_lf / total_num),
        x = 0.84, y = 0.42, size = 28,
        colour='#33CC00', alpha=0.75
    ) +
    draw_text(
        sprintf('%.1f%%', 100 * passed_nf / total_num),
        x = 0.84, y = 0.3, size = 28,
        colour='#33CC00', alpha=0.75
    ) +
    draw_text(
        sprintf('Total: %.1f%%', 100 * passed_all / total_num),
        x = 0.84, y = 0.05, size = 28,
        colour='#33CC00', alpha=0.75
    )
# print(p1)
print(p2)
dev.off()

# stat.s
total           <- nrow(meta)
total_w_asm     <- length(unique(df_melt$ID))
total_rel_taxa  <- length(unique(df_melt$ID[df_melt$Taxon != 'Other']))
total_pgs       <- length(unique(df_melt$ID[df_melt$InPG == 'Yes']))
write(sprintf(
    'Plot over all samples with assembly: %d samples (%.1f of all)',
    total_w_asm, 100 * total_w_asm / total
), file=log_file, append=TRUE)
write(sprintf(
    'Samples from relevant taxa (w/o \"Other\"): %d samples (%.1f of all, %.1f of with assembly)',
    total_rel_taxa, 100 * total_rel_taxa / total, 100 * total_rel_taxa / total_w_asm
), file=log_file, append=TRUE)
write(sprintf(
    'Plot over all samples from pan-genomes: %d samples (%.1f of all, %.1f of with assembly)',
    total_pgs, 100 * total_pgs / total, 100 * total_pgs / total_w_asm
), file=log_file, append=TRUE)
# Passing assembly quality criteria
passed_all <- sum(df$AssemblyNum <= CUTOFF_CONTIG & df$AssemblyN50 >= CUTOFF_N50 & df$AssemblyL50 <= CUTOFF_L50, na.rm=TRUE)
passed_pg  <- sum(df$AssemblyNum <= CUTOFF_CONTIG & df$AssemblyN50 >= CUTOFF_N50 & df$AssemblyL50 <= CUTOFF_L50 & df$InPG == "Yes", na.rm=TRUE)
write(sprintf('Passing quality criteria: %d (%.1f of %d w/ assembly); %d in pan-genomes', passed_all, 100 * passed_all / nrow(df), nrow(df), passed_pg), log_file, append=TRUE)
# GC%
write('\nMedian GC content: all vs. in pan-genome', log_file, append=TRUE)
for(taxon in sort(setdiff(unique(df_melt$Taxon),'Other'))){
    gc_all <- median(df_melt[df_melt$Taxon==taxon & df_melt$variable=='AssemblyGC','value'])
    gc_pg  <- median(df_melt[df_melt$Taxon==taxon & df_melt$variable=='AssemblyGC' & df_melt$InPG == 'Yes','value'])
    write(sprintf('\t%s: %.1f / %.1f', taxon, gc_all, gc_pg), log_file, append=TRUE)
}
