#!/usr/bin/Rscript

## rpart.plot: prp function:
##  http://www.milbo.org/rpart-plot/
##  http://www.milbo.org/rpart-plot/prp.pdf
##  https://cran.r-project.org/web/packages/rpart.plot/rpart.plot.pdf

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(rpart)) # trees
suppressMessages(library(rpart.plot)) # for plotting

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    # required parameters
    parser$add_argument('--tree_rds',   '-t', help='', required=TRUE, nargs='+')
    parser$add_argument('--tree_annot', '-a', help='', required=TRUE, nargs='+')
    parser$add_argument('--tree_name',  '-n', help='', required=TRUE, nargs='+')
    parser$add_argument('--ofile',      '-o',  help='', required=TRUE)
    # prp args
    parser$add_argument('--prp_type',   default=0, type='integer')
    parser$add_argument('--prp_extra',  default=0, type='integer')
    parser$add_argument('--prp_yesno',  default=1, type='integer')
    parser$add_argument('--prp_left',   default=TRUE)
    parser$add_argument('--prp_varlen', default=0, type='integer')
    parser$add_argument('--prp_faclen', default=0, type='integer')
    parser$add_argument('--prp_cex')
#     parser$add_argument('--prp_', default=)
    # pdf args
    parser$add_argument('--pdf_width',  default=9, type='integer')
    parser$add_argument('--pdf_height', default=6, type='integer')
    # other
    parser$add_argument('--verbose', '-v', action='store_true')
    return(parser)
}

# Plot tree
## TODO
draw_tree <- function(x, type, extra, yesno, left, varlen, faclen, cex, node.fun, split.fun, box.palette, main){
    tryCatch({
        do.call('prp', as.list(match.call()[-1]))
    },
        warning=function(w) { write( sprintf('WARNING:\n%s\nfor \"%s\"',w,main), stdout()) },
        error=function(e)   { plot.new(); write(sprintf('ERROR:\n%s\nfor \"%s\"',e,main,sep=''), stdout()) },
        finally={}
    )
}

# custom node labels: 0/1 -> susceptible/resistant
node_bin2str <- function(x, labs, digits, varlen){
    labs <- sub('0\n','susceptible\n',labs)
    labs <- sub('1\n','resistant\n',labs)
    return(labs)
}

# custom split labels: 0/1 -> absent/present
split_genes_bin2str <- function(x, labs, digits, varlen, faclen){
    labs <- sub(' = 0$',' is absent\n',labs)
    labs <- sub(' = 1$',' is present\n',labs)
    return(labs)
}

# custom split labels: 0/1 -> alt/ref
split_snps_bin2str <- function(x, labs, digits, varlen, faclen){
    labs <- sub(' = 0$',' is ref. allele\n',labs)
    labs <- sub(' = 1$',' is alt. allele\n',labs)
    return(labs)
}

# custom split labels: feature ID -> feature annotation
split_annot <- function(x, labs, digits, varlen, faclen, annot){
    labs <- sapply(labs, function(s){
        s <- unlist(strsplit(s,' '))
        s[1] <- annot[s[1]]
        return(paste(s, collapse=' '))
    })
    return(labs)
}


## MAIN
## Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))
testit::assert(length(args$tree_annot)==length(args$tree_rds))
testit::assert(length(args$tree_annot)==length(args$tree_name))

## Plot
pdf(args$ofile, width=args$pdf_width, height=args$pdf_height)
for(i in 1:length(args$tree_rds)){
    # tree files: model and annotation
    f_t <- args$tree_rds[i]
    f_a <- args$tree_annot[i]
    # name, model, annotation
    i_n <- args$tree_name[i]
    i_t <- readRDS(f_t)
    i_a <- read.csv(file=f_a, header=TRUE, sep='\t', check.names=FALSE, stringsAsFactors=FALSE)
    rownames(i_a) <- i_a$Feature
    
    # get relevant annotation information
    if('ref' %in% colnames(i_a)){ # variant
        i_a_ <- apply(i_a[,c('product','DNA_change')], 1, paste, collapse='\n')
        names(i_a_) <- rownames(i_a)
        split_fun <- function(x, labs, digits, varlen, faclen) split_annot(x, split_snps_bin2str(x, labs, digits, varlen, faclen), digits, varlen, faclen, annot=i_a_)
    } else { # gene/centroid
        i_a_ <- i_a$product
        names(i_a_) <- rownames(i_a)
        split_fun <- function(x, labs, digits, varlen, faclen) split_annot(x, split_genes_bin2str(x, labs, digits, varlen, faclen), digits, varlen, faclen, annot=i_a_)
    }
    
    draw_tree(
        x=i_t,
        type=args$prp_type, extra=args$prp_extra, yesno=args$prp_yesno,
        left=args$prp_left, varlen=args$prp_varlen, faclen=args$prp_faclen, cex=args$prp_cex,
        node.fun=node_bin2str,
        split.fun=split_fun,
        box.palette="-Greens",
        main=i_n
    )
}
dev.off()
