## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(ape))
# suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(RColorBrewer))

# /share/runs/Roary/no_split/pangenome_analysis
# args <- list(
#     tree='essential_genes_core_raxml/RAxML_result.egenes_core_tree',
#     meta_rds=c('sample_data.rds', 'sample_data_rest.rds'),
#     src_path="~/git_repos/Bacteria/SubProject_Pangenomes/scripts"
# )

## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--tree', help='', required=TRUE)
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
# Tree
raxml <- ape::read.tree(args$tree) # read in
raxml <- ape::makeNodeLabel(raxml) # add node labels
write(sprintf('Read in tree with %d tips\n', length(raxml$tip.label)), log_file, append=TRUE)

# Sample meta-data
tmp <- read_meta_rds(args$meta_rds)
meta <- tmp$meta
samples <- tmp$samples

# Tree meta data
raxml_meta <- meta[raxml$tip.label,] # subset of samples
raxml_meta$Color <- convert_taxa(raxml_meta$TaxSpecies) # color = taxon
raxml_meta$Label <- abbr_taxa(convert_taxa(raxml_meta$TaxSpecies)) # label = abbr. taxon
# remove some labels
last_org <- NULL
for(tip in raxml$tip.label){
    if(is.null(last_org)){
        last_org <- list(Label=raxml_meta[tip, 'Label'], Tips=c(tip))
        next
    }
    if(last_org$Label == raxml_meta[tip, 'Label']){
        last_org$Tips <- c(last_org$Tips, tip)
    } else {
        # print(sprintf('Old: %s (%d) | New: %s', last_org$Label, length(last_org$Tips), raxml_meta[tip, 'Label']))
        write(sprintf('Old: %s (%d) | New: %s', last_org$Label, length(last_org$Tips), raxml_meta[tip, 'Label']), log_file, append=TRUE)
        # Set label for the tip in the "middle", other to empty string
        raxml_meta[last_org$Tips,'Label'] <- ''
        raxml_meta[last_org$Tips[round(length(last_org$Tips)/2)],'Label'] <- last_org$Label
        # new organism
        last_org <- list(Label=raxml_meta[tip, 'Label'], Tips=c(tip))
    }
    # last tip
    if(tip == raxml$tip.label[length(raxml$tip.label)]){
        # print(sprintf('Last: %s (%d)', last_org$Label, length(last_org$Tips)))
        write(sprintf('Last: %s (%d)', last_org$Label, length(last_org$Tips)), log_file, append=TRUE)
        if(last_org$Label %in% raxml_meta$Label){
            # Remove same label from earlier found organism
            raxml_meta[raxml_meta$Label == last_org$Label,'Label'] <- ''
        }
        # Set label for the tip in the "middle", other to empty string
        raxml_meta[last_org$Tips,'Label'] <- ''
        raxml_meta[last_org$Tips[round(length(last_org$Tips)/2)],'Label'] <- last_org$Label
    }
}
# Add leading whitespace in labels (because offset does not work)
raxml_meta$Label <- paste(paste(rep(' ', 2),collapse=''), raxml_meta$Label, sep='')
# Cast color labels to ordered factors -> better color assignment
raxml_meta$Color <- factor(x=raxml_meta$Color, levels=sort(unique(raxml_meta$Color)), ordered=TRUE)

## Plot
# TODO: ggtree seems to ignore offset in geom_tiplab

pdf(sprintf('%s.pdf', args$obname), height=6, width=6)

# Plot whole tree
p <- ggtree(raxml, layout='circular')
p <- p %<+% raxml_meta +
    geom_tiplab(size=2, offset=1, aes(angle=angle, label=Label)) +
    geom_tippoint(size=1, shape=20, alpha=0.5, aes(color=Color)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(18))
print(p)

# p <- ggtree(raxml, layout='rectangular')
# p <- p %<+% raxml_meta +
#     geom_tiplab(size=3, offset=1, aes(label=Label)) +
#     geom_tippoint(size=0.5, shape=15, alpha=0.5, aes(color=Color)) +
#     scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(18))
# print(p)

# # Plot sub-tree(s) per taxon
# for(taxon in sort(unique(raxml_meta$Color))){
#     # gzoom(raxml, grep(taxon, raxml_meta$Color))
#     # tips of that taxon
#     taxon_tips <- raxml$tip.label[ grepl(taxon, raxml_meta$Color) ]
#     # LCA
#     taxon_mrca <- getMRCA(raxml, taxon_tips)
#     # Complete tree
#     p <- ggtree(raxml, layout='circular') + geom_tippoint() + ggtitle(taxon)
#     # Get position of LCA
#     cpos <- get_clade_position(p, node=taxon_mrca)
#     # Get subtree
#     p_taxon <- with(cpos, p + xlim(xmin, xmax) + ylim(ymin, ymax))
#     print(p_taxon)
# }
dev.off()
