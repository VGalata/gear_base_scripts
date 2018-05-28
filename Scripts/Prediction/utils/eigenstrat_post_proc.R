#!/usr/bin/Rscript

## Process binary feature file
## Assumption: values are 0/1, missing NA (but can be specified)

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(parallel))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('-eig_chisq', help='EIGENSTRAT Chisq results', required=TRUE)
    parser$add_argument('-eig_snps', help='SNP file used as input for EIGENSTRAT', required=TRUE)
    parser$add_argument('-eig_ev', help='EIGENSTRAT PCA Eigenvalues', required=TRUE)
    parser$add_argument('-eig_k', help='Used number of PCs in EIGENSTRAT PCA', required=TRUE, type='integer')
    # if performance computation should be performed
    parser$add_argument('-geno_file', help='Optional')
    parser$add_argument('-pheno_file', help='Optional')
    parser$add_argument('-pheno_name', help='Optional')
    parser$add_argument('-samples', help='Optional')
    parser$add_argument('-features', help='Optional')
    #
    parser$add_argument('--adj', default='fdr', choices=p.adjust.methods)
    parser$add_argument('--alpha', default=0.01, type='double')
    parser$add_argument('--sel_ofile', default='')
    #
    parser$add_argument('--src_path', help='', default='.')
    parser$add_argument('--cores', help='', default=2, type='integer')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    
    return(parser)
}

## MAIN
# Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
# Print info
if(args$verbose){
    write(sprintf('\nArguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

## LIBS
src_path <- args$src_path
source(sprintf("%s/my_utils.R",src_path))
source(sprintf("%s/perf_utils.R",src_path))
source(sprintf("%s/pv_qq_plot.R",src_path))
source(sprintf("%s/pc_plot.R",src_path))

## READ IN
# association statistics
eig_chisq <- read.csv(file=args$eig_chisq, header=TRUE, sep=' ')
# used features
eig_snps  <- read.csv(file=args$eig_snps, header=FALSE, sep='\t', stringsAsFactors=FALSE)
eig_snps  <- unlist(eig_snps[,1])
eig_snps  <- gsub('^rs','',eig_snps)
# rownames = features
rownames(eig_chisq) <- eig_snps
# add feature as column
eig_chisq$Feature   <- eig_snps

if(args$verbose){
    write(sprintf('Read in file has %d columns and %d rows',ncol(eig_chisq),nrow(eig_chisq)),stdout())
    print(head(eig_chisq))
}

## P-values and adjustment
if(! 'Chisq_pvalue' %in% colnames(eig_chisq)){
    eig_chisq$Chisq_pvalue     <- pchisq(eig_chisq$Chisq, df=1, lower.tail=FALSE, log.p=FALSE)
}
if(! 'EChisq_pvalue' %in% colnames(eig_chisq)){
    eig_chisq$EChisq_pvalue    <- pchisq(eig_chisq$EIGENSTRAT, df=1, lower.tail=FALSE, log.p=FALSE) # see EIGENSTRAT paper df = 1
}
if(! 'Chisq_pvalue_adj' %in% colnames(eig_chisq)){
    eig_chisq$Chisq_pvalue_adj <- p.adjust(p=eig_chisq$Chisq_pvalue, method=args$adj)
}
if(! 'EChisq_pvalue_adj' %in% colnames(eig_chisq)){
    eig_chisq$EChisq_pvalue_adj<- p.adjust(p=eig_chisq$EChisq_pvalue, method=args$adj)
}

## Plots
## Q-Q-Plot
pdf(sprintf("%s.pdf",args$eig_chisq), width=5, height=5)
ggd.qqplot(eig_chisq$Chisq_pvalue,      sprintf("QQ-plot of Chi-sq. p-values"))
ggd.qqplot(eig_chisq$Chisq_pvalue_adj,  sprintf("QQ-plot of Chi-sq. p-values (adj.)"))
ggd.qqplot(eig_chisq$EChisq_pvalue,     sprintf("QQ-plot of EIGENSTRAT p-values"))
ggd.qqplot(eig_chisq$EChisq_pvalue_adj, sprintf("QQ-plot of EIGENSTRAT p-values (adj.)"))
## PC var. http://stats.stackexchange.com/questions/31908/what-is-percentage-of-variance-in-pca
ev    <- unlist( read.csv(file=args$eig_ev, header=FALSE, stringsAsFactors=FALSE) ) # should remove leading whitespaces and cast to numeric
p <- plot_pc_var(ev=ev, k=args$eig_k)
print(p)
dev.off()

## Sort by adj. p-value of adj. chisq.
eig_chisq <- eig_chisq[with(eig_chisq,order(EChisq_pvalue_adj)),]

## Add features from same cluster
if(!is.null(args$eig_snps_cl)){
    snps_cl <- read.csv(file=args$eig_snps_cl, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    eig_chisq$SameCluster <- sapply(rownames(eig_chisq), function(snp){
        snp_cl <- snps_cl$Cluster[snps_cl$ID==snp]
        return( paste(sort(setdiff(snps_cl$ID[snps_cl$Cluster==snp_cl],snp)),collapse=';') )
    })
}

## Compute perf. if geno/pheno matrices given
if(!is.null(args$geno_file) & !is.null(args$pheno_file) & !is.null(args$pheno_name)){
    # read in
    x <- read.csv(file=args$geno_file, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
    x <- t(x)
    y <- read.csv(file=args$pheno_file, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
    # select samples
    if(!is.null(args$samples)){
        include_samples <- read_list(list_file=args$samples, verbose=args$verbose)
        x <- sel_mat_rows(mat=x, incl_rnames=include_samples, verbose=args$verbose)
        y <- sel_mat_rows(mat=y, incl_rnames=include_samples, verbose=args$verbose)
        # all samples to include should be there
        testit::assert(sprintf(
            "Not all samples which should be included are there:\n%s\n",
            paste(include_samples[ sapply(include_samples, function(n){ !(n %in% rownames(x)) }) ], collapse=',')
        ), all( sapply(include_samples, function(n){ n %in% rownames(x) }) ))
        testit::assert(sprintf(
            "Not all samples which should be included are there:\n%s\n",
            paste(include_samples[ sapply(include_samples, function(n){ !(n %in% rownames(y)) }) ], collapse=',')
        ), all( sapply(include_samples, function(n){ n %in% rownames(y) }) ))
    }
    # select features
    if(!is.null(args$features)){
        include_features <- read_list(list_file=args$features, verbose=args$verbose)
        x <- sel_mat_cols(mat=x, incl_cnames=include_features, verbose=args$verbose)
    }
    # select pheno
    y_old <- y
    y <- y[,args$pheno_name,drop=TRUE]
    names(y) <- rownames(y_old)
    # sort by sample ID
    x <- x[sort(rownames(x)),]
    y <- y[sort(names(y))]
    # check
    testit::assert(all(rownames(x)==names(y)))

    # perf
    perf_num <- c('N','NX','NY')
    perf_flag <- 'flag'
    perf_corr <- 'cor'
    perf_conf <- c('TN','TP','FN','FP')
    perf_perf <- c('ERR','ACC','B_ACC','SENS','SPEC','PREC','NPV','FPR','FNR','Fmeasure','gm_RS','gm_RP')
    perf_m <- c(perf_num, perf_corr, perf_flag, perf_conf, perf_perf)
    
    eig_chisq <- data.frame(
        eig_chisq,
        matrix(NA, nrow=nrow(eig_chisq), ncol=length(perf_m), dimnames=list(rownames(eig_chisq),perf_m)),
        row.names=rownames(eig_chisq), stringsAsFactors=FALSE, check.names=FALSE
    )
    
    successed <- FALSE; tried <- 0
    while(!successed & tried <5){
        cl <- tryCatch(makeCluster(args$cores, type = "FORK"), error=function(e){NULL} )
        if(is.null(cl)){ tried <- tried + 1; next }
        
#         clusterExport(cl=cl, varlist=c('x','y','src_path'), envir=.GlobalEnv)
        
        eig_chisq$N     <- rep(length(y), nrow(eig_chisq))
        eig_chisq$NX    <- parSapply(cl=cl, X=rownames(eig_chisq), FUN=function(snp){ paste(sum(x[,snp]==0, na.rm=TRUE),sum(x[,snp]==1, na.rm=TRUE), sep=':') }, USE.NAMES=TRUE)
        eig_chisq$NY    <- rep(paste(sum(y==0, na.rm=TRUE),sum(y==1, na.rm=TRUE),sep=':'), nrow(eig_chisq))
        eig_chisq$cor   <- parSapply(cl=cl, X=rownames(eig_chisq), FUN=function(snp){ cor(x=na.omit(x[,snp]), y=y[!is.na(x[,snp])], method='pearson') }, USE.NAMES=TRUE)
        eig_chisq$flag  <- sapply(eig_chisq$cor, function(snp_cor){ ifelse(!is.na(snp_cor) & snp_cor<0,'abs->res','prs->res') })
        
        res_conf <- parSapply(cl=cl, X=rownames(eig_chisq), FUN=function(snp){
            source(sprintf("%s/perf_utils.R",src_path))
            snp_pred <- x[,snp]
            snp_cor <- cor(x=na.omit(snp_pred), y=y[!is.na(snp_pred)], method='pearson')
            if( !is.na(snp_cor) & (snp_cor < 0) ){ # neg. corr -> reverse pred.
                snp_pred <- sapply(snp_pred, function(p){ ifelse(p==1,0,ifelse(p==0,1,p)) })
            }
            cm <- confusion_mat(true=y, pred=snp_pred, levels=c(0,1))
            return(c(TN=cm$TN,TP=cm$TP,FN=cm$FN,FP=cm$FP))
        }, USE.NAMES=TRUE)
        eig_chisq[,perf_conf] <- t(res_conf)[rownames(eig_chisq),perf_conf]
        
        res_perf <- parSapply(cl=cl, X=rownames(eig_chisq), FUN=function(snp){
            source(sprintf("%s/perf_utils.R",src_path))
            snp_pred <- x[,snp]
            snp_cor <- cor(x=na.omit(snp_pred), y=y[!is.na(snp_pred)], method='pearson') #cor(x=snp_pred, y=y, method='pearson')
            if( !is.na(snp_cor) & (snp_cor < 0) ){ # neg. corr -> reverse pred.
                snp_pred <- sapply(snp_pred, function(p){ ifelse(p==1,0,ifelse(p==0,1,p)) })
            }
            return( all_perf(true=y, pred=snp_pred, levels=c(0,1), with_names=TRUE) )
        }, USE.NAMES=TRUE)
        eig_chisq[,perf_perf] <- t(res_perf)[rownames(eig_chisq),perf_perf]
        
        stopCluster(cl)
        successed <- TRUE
    }
    if(!successed){ stop(sprintf("Could not perform parallel computing")) }
}

## Save
write.table(x=eig_chisq, file=sprintf("%s.adj",args$eig_chisq), row.names=FALSE, sep='\t', quote=FALSE)

## If alpha given: list of features passed the threshold
if(!is.null(args$alpha)){
    sel_snps <- eig_chisq$Feature[!is.na(eig_chisq$EChisq_pvalue_adj) & eig_chisq$EChisq_pvalue_adj <= args$alpha]
    if(args$verbose){
        write(sprintf(
            "From %d tested features: %d with NA, %d selected",
            nrow(eig_chisq), sum(is.na(eig_chisq$EChisq_pvalue_adj)), length(sel_snps)
        ),stdout())
    }
    if(args$sel_ofile==""){ args$sel_ofile <- sprintf("%s.sel",args$eig_chisq) }
    write(x=sel_snps, file=args$sel_ofile)
}