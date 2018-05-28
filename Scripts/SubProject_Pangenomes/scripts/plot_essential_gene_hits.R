#!/usr/bin/Rscript

mcn_cutoff <- 1.25

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
    parser$add_argument('--csmap', help='', required=TRUE)
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
# Centroid/subject mapping
csmap <- read.csv(file=args$csmap, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
# essential gene annotations
annot <- read.csv(file=args$annot, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
rownames(annot) <- annot$Name
# sample meta data
tmp <- read_meta_rds(args$meta_rds)
meta <- tmp$meta
samples <- tmp$samples
write(sprintf('Meta data: %d samples (%d, rest %d)', nrow(meta), length(samples$main), length(samples$rest)), log_file, append=TRUE)

# Aggregate subject based mapping: per subject = mean #centroids per pan-genome
csmap$SampleID <- sapply(csmap$CentroidID, geneID2sampleID) # get sample ID from centroid ID
csmap <- csmap[sapply(csmap$SampleID, function(x){ x %in% samples$main }),] # keep only relevant samples
testit::assert(length(setdiff(unique(csmap$SampleID), samples$main)) == 0) # check
csmap$Taxon <- convert_taxa(meta[csmap$SampleID,'TaxSpecies']) # get sample taxon
csmap_aggr <- aggregate(csmap$CentroidID, by=list(SubjectID=csmap$SubjectID, Taxon=csmap$Taxon), FUN=length) # num. centroids per taxon/pan-genome)
csmap_aggr_pg <- csmap_aggr # backup
csmap_aggr <- aggregate(csmap_aggr$x, by=list(SubjectID=csmap_aggr$SubjectID), FUN=mean) # mean num. centroids over pan-genomes per subject
rownames(csmap_aggr) <- csmap_aggr$SubjectID

# Log stat.s
csmap_exceptions <- data.frame(
    EssentialGeneID=csmap_aggr[csmap_aggr$x >= mcn_cutoff,'SubjectID'],
    EssentialGeneAccession=annot[csmap_aggr[csmap_aggr$x >= mcn_cutoff,'SubjectID'],'Accession'],
    EssentialGeneDescription=annot[csmap_aggr[csmap_aggr$x >= mcn_cutoff,'SubjectID'],'Description'],
    MeanCentroidNumber=csmap_aggr[csmap_aggr$x >= mcn_cutoff,'x'],
    stringsAsFactors=FALSE
)
write.table(x=csmap_exceptions, file=sprintf('%s_mult_centroids.tsv', args$obname), sep='\t', row.names=FALSE, quote=FALSE)
write(sprintf('Essential genes with mean centroid number >= %f', mcn_cutoff), log_file, append=TRUE)
for(i in 1:nrow(csmap_aggr)){
    mean_cnum <- csmap_aggr[i,'x']
    s_id <- csmap_aggr[i,'SubjectID']
    if(mean_cnum >= mcn_cutoff){
        write(sprintf('%s (%s): %f', s_id, annot[s_id,'Description'], mean_cnum), log_file, append=TRUE)
        tmp <- csmap_aggr_pg[csmap_aggr_pg$SubjectID == s_id,]
        for(j in 1:nrow(tmp)){
            j_taxon <- tmp[j,'Taxon']
            write(sprintf('\t%s: %d: %s', j_taxon, tmp[j,'x'], paste(csmap[csmap$Taxon==j_taxon & csmap$SubjectID==s_id,'CentroidID'], collapse=',')), log_file, append=TRUE)
        }
    }
}

# In XX% of all samples
XXpct <- 99
in_XXpct <- sum(objs$ov_g_s_prs >= XXpct)
write(sprintf('Essential genes found in >=%f pct of isolates: %d (%.1f of %d)', XXpct, in_XXpct, 100 * in_XXpct / length(objs$ov_g_s_prs), length(objs$ov_g_s_prs)), log_file, append=TRUE)

# Re-shape: taxon, essential gene, count -> taxon x essential gene
csmap_aggr_pg <- dcast(csmap_aggr_pg, Taxon ~ SubjectID, drop=FALSE, fill=0) # re-shape
rownames(csmap_aggr_pg) <- csmap_aggr_pg$Taxon # set row names
csmap_aggr_pg <- csmap_aggr_pg[,setdiff(colnames(csmap_aggr_pg), 'Taxon')] # rm column with row names
testit::assert(all(sort(rownames(csmap_aggr_pg)) ==  rownames(objs$ov_g_sp)))
for(s_id in colnames(objs$ov_g_sp)){ # add missing essential genes (those without any hits)
    if(! s_id %in% colnames(csmap_aggr_pg)){
        csmap_aggr_pg[,s_id] <- 0
    }
}

## NOTE Plots
write('Create plots', log_file, append=TRUE)
pdf(sprintf('%s.pdf', args$obname), width=9, height=6)
## Plot: Samples x essential genes
# superheat(
#     X=objs$ov_g_s,
#     pretty.order.rows=TRUE,
#     pretty.order.cols=TRUE,
#     row.dendrogram=FALSE,
#     col.dendrogram=FALSE,
#     smooth.heat=FALSE,
#     scale=FALSE,
#     yt = round(objs$ov_g_s_prs),
#     yt.axis.name = "Presence [%]",
#     yt.plot.type = "bar",
#     yt.obs.col = sapply(round(objs$ov_g_s_prs), function(x){ ifelse(x<objs$min_prs,'#CCCCCC','#666666') }),
#     heat.pal=c("#3B528BFF", "#5DC863FF", "#FDE725FF"), # subset of viridis scheme colors
#     heat.lim = c(1, max(objs$ov_g_s)), # 0 = missing (will get NA color)
#     heat.na.col = "#CCCCCC",
#     legend=TRUE,
#     grid.hline.size=0.1, grid.vline.size=0.1,
#     # bottom.label='variable', bottom.label.text.angle='90', bottom.label.text.size=3,
#     # left.label.text.size=3, left.label.size=0.5,
#     row.title='Samples', column.title='Essential genes'
# )

ov_g_s_ss <- names(objs$ov_g_s_prs)[objs$ov_g_s_prs < objs$min_prs] # only those with prs. less than min cutoff
ov_g_s_taxa <- abbr_taxa(convert_taxa(meta[rownames(objs$ov_g_s),'TaxSpecies']))
u_taxa <- sort(unique(ov_g_s_taxa))
s <- paste(rep(' ', 25), collapse='')
ov_g_s_taxa <- sapply(ov_g_s_taxa, function(x){ ifelse((which(u_taxa==x) %% 2) == 0, paste(x, s, sep=''), paste(s, x, sep='')) })
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
    yt.obs.col = rep('#666666', length(ov_g_s_ss)), # sapply(objs$ov_g_s_prs[ov_g_s_ss], function(x){ ifelse(x>=95, '#FF3333', '#333333') }),
    # heat.pal=c("#3B528BFF", "#5DC863FF", "#FDE725FF"), # subset of viridis scheme colors
    # heat.lim = c(1, max(objs$ov_g_s[,ov_g_s_ss])), # 0 = missing (will get NA color)
    # heat.na.col = "#CCCCCC",
    heat.col.scheme='viridis',
    heat.pal.values=my_pal_values,
    legend.breaks=c(0, 1, 2, 5),
    legend=TRUE,
    grid.hline.size=0.1, grid.vline.size=0.1,
    bottom.label='variable', bottom.label.text.angle='90', bottom.label.text.size=3, bottom.label.size=0.35,
    left.label.text.size=3, left.label.size=0.3,
    row.title='Samples', column.title='Essential genes', legend.title='# hits'
)

## Plot: Per species and cat. boxplot of total/unique/mult. hits
pt <-
    theme(
        axis.ticks=element_blank(),
        panel.grid.major=element_line(colour="#CCCCCC"),
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

custom_labeller <- c(
    total="Total hits",
    uniq="Unique hits",
    mult="Multiple hits"
)

p <-ggplot(data=objs$ov_n_s, aes(x=taxon, y=value)) +
    geom_boxplot() +
    facet_wrap(~variable, ncol=1, scales='free_y', labeller=labeller(variable=custom_labeller)) +
    xlab('') + ylab('Number of hits') +
    pt
print(p)

## Plot: Species x essential genes

# egene_group <- function(egene_name){
#     if(grepl('(R|r)ibosomal protein', egene_name)){
#         return('Ribosomal protein')
#     }
#     if(grepl('tRNA ligase', egene_name)){
#         return('tRNA ligase')
#     }
# }

# species x genes = mean #hits
my_pal_values <- scales::rescale(c(0,0.5,0.75,1,max(objs$ov_g_sp)), to=c(0,1))
superheat(
    X=objs$ov_g_sp,
    pretty.order.rows=TRUE,
    pretty.order.cols=TRUE,
    row.dendrogram=TRUE,
    col.dendrogram=FALSE,
    linkage.method="average",
    smooth.heat=FALSE,
    scale=FALSE,
    legend=TRUE,
    yt = round(objs$ov_g_s_prs[colnames(objs$ov_g_sp)]),
    yt.axis.name = "Presence [%]",
    yt.plot.type = "bar",
    # yt.obs.col = sapply(objs$ov_g_s_prs[colnames(objs$ov_g_sp)], function(x){ ifelse(x>=99, '#0099cc', '#333333') }),
    heat.col.scheme='viridis',
    heat.pal.values=my_pal_values,
    legend.breaks=c(0, 1, 2, 5),
    grid.hline.size=0.1, grid.vline.size=0.1,
    left.label.text.size=3, left.label.size=0.5,
    row.title='Taxa', column.title='Essential genes', legend.title='Mean # hits'
)

# species x genes = mean #centroids
my_pal_values <- scales::rescale(c(0,0.5,0.75,1,max(csmap_aggr_pg)), to=c(0,1))
superheat(
    X=csmap_aggr_pg,
    pretty.order.rows=TRUE,
    pretty.order.cols=TRUE,
    row.dendrogram=FALSE,
    col.dendrogram=FALSE,
    linkage.method="average",
    smooth.heat=FALSE,
    scale=FALSE,
    legend=TRUE,
    yt = csmap_aggr[colnames(csmap_aggr_pg),'x'],
    yt.obs.col = sapply(csmap_aggr[colnames(csmap_aggr_pg),'x'],function(x){ ifelse(x >= mcn_cutoff, '#FF3333', '#333333') }),
    yt.axis.name = "Mean number\nof centroid hits",
    yt.plot.type = "bar",
    heat.col.scheme='viridis',
    heat.pal.values=my_pal_values,
    legend.breaks=c(0, 1, 2, 5),
    grid.hline.size=0.1, grid.vline.size=0.1,
    left.label.text.size=3, left.label.size=0.5,
    row.title='Taxa', column.title='Essential genes', legend.title='Mean # centroids'
)

# superheat(
#     X=objs$ov_g_sp[,apply(objs$ov_g_sp, 2, function(x){ any(x < 0.5) })],
#     pretty.order.rows=TRUE,
#     pretty.order.cols=TRUE,
#     row.dendrogram=TRUE,
#     col.dendrogram=FALSE,
#     smooth.heat=FALSE,
#     scale=FALSE,
#     legend=TRUE,
#     grid.hline.size=0.1, grid.vline.size=0.1,
#     bottom.label='variable', bottom.label.text.angle='90', bottom.label.text.size=3, bottom.label.size=0.4,
#     left.label.text.size=3, left.label.size=0.5,
#     row.title='Taxa', column.title='Essential genes'
# )

dev.off()
