#!/usr/bin/Rscript

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
# suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
# suppressMessages(library(splines))
# suppressMessages(library(MASS))
suppressMessages(library(gtools))

# test args
# args <- list(
#     tabs=list.files(path='..', pattern='pan_genome_perm.tsv', full.names=TRUE, recursive=TRUE),
#     src_path="~/git_repos/Bacteria/SubProject_Pangenomes/scripts"
# )

## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--tabs', help='', required=TRUE, nargs='+')
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
read_in_perms <- function(tab_file){
    tab <- read.csv(file=tab_file, header=TRUE, sep='\t', check.names=FALSE, comment.char = "#", stringsAsFactors=FALSE)
    # Remove core of 100% (not informative - goes down with increasing number of genomes)
    tab <- tab[,colnames(tab)[ !grepl('NumCore100',colnames(tab)) ]]
    # Change column names
    colnames(tab) <- sapply(colnames(tab), function(n){
        if(n=='NumAll'){  return('Total') }
        if(n=='NumNew'){  return('New') }
        if(n=='NumUniq'){ return('Unique') }
        if(grepl('NumCore',n)){
            return( paste(gsub('Core','Core ',gsub('Num','',n)),'%',sep='') )
        } else {
            return(n)
        }
    })
    # Re-shape
    tab <- reshape::melt.data.frame(tab, id.vars=c('Perm','N'))
    # Aggregate: compute medians
    tab <- aggregate(tab$value, by=list(N=tab$N, variable=tab$variable), FUN=median)
    return(tab)
}

df <- NULL
for(tab_file in args$tabs){
    tab <- read_in_perms(tab_file)
    tab_taxon <- gsub('_', ' ', basename(dirname(dirname(tab_file))))
    tab$Taxon <- rep(tab_taxon, nrow(tab))
    if(is.null(df)){
        df <- tab
    } else {
        df <- rbind(df, tab)
    }
}

## Taxa
taxa <- sort(unique(df$Taxon))

## Heaps' Law
hl <- data.frame(
    Taxon=taxa,
    eq=NA, a=NA, b=NA,
    stringsAsFactors=FALSE, check.names=FALSE, row.names=taxa
)
for(taxon in taxa){
    total_nls <- nls(x ~ a*N^b, start=list(a=1,b=1), data=df[df$Taxon==taxon & df$variable=='Total',], control=nls.control(minFactor=1e-10))
    hl[taxon,'eq']<- sprintf('n = %.2f * N ^ %.2f',coef(total_nls)[1],coef(total_nls)[2])
    hl[taxon,'a'] <- coef(total_nls)[1]
    hl[taxon,'b'] <- coef(total_nls)[2]
}
write.table(x=hl, file=sprintf('%s.tsv', args$obname), row.names=FALSE, sep='\t', quote=FALSE)

## Plot data
linetypes   <- c('solid', 'dashed', 'dotted')
linecolours <- scales::hue_pal()(6)
line_params <- expand.grid(linecolours, linetypes)
df_lt <- sapply(taxa, function(x){ as.character(line_params[which(taxa == x),2]) })
df_lc <- sapply(taxa, function(x){ as.character(line_params[which(taxa == x),1]) })


## Plot
pt <-
    theme(
        axis.ticks = element_blank(),
        panel.grid.major=element_line(colour="#CCCCCC"),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=16, margin = margin(t = 10, b = 10)),
        axis.title.y = element_text(size=16, margin = margin(l = 10, r = 10)),
        axis.text=element_text(colour="black"),
        panel.background=element_rect(fill="white"),
        panel.border = element_rect(colour="#333333", fill=NA),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.key.width = unit(2.5,'line'),
        legend.text = element_text(size=10)
    )

pdf(sprintf('%s.pdf', args$obname), width=9, height=6)
for(var_type in levels(df$variable)){
    p <-ggplot(subset(df, variable == var_type), aes(x=N, y=x)) +
        scale_linetype_manual(values=df_lt, name='') +
        scale_colour_manual(values=df_lc, name='') +
        geom_line(size=1, aes(linetype=Taxon, colour=Taxon)) +
        scale_x_sqrt() +
        ylab("Number of centroids") +
        xlab("Number of genomes") +
        ggtitle(var_type) +
        pt
    if(var_type == 'New'){
        p <- p + scale_y_sqrt()
    } else {
        p <- p +
        scale_y_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        )
    }
    print(p)
}
dev.off()
