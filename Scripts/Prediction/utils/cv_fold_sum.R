#!/usr/bin/Rscript

## LIBS
suppressMessages(library(argparse))
suppressMessages(library(testit))
suppressMessages(library(tools))



## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    # required parameters
    parser$add_argument('-pred_files', help='', required=TRUE, nargs="+")
    parser$add_argument('-o_file', help='Output file', required=TRUE)
    # other
    parser$add_argument('--plot', help='', action='store_true')
    parser$add_argument('--src_path',   help='Print additional information.', default='.')
    parser$add_argument('--verbose', help='Print additional information.', action='store_true')
    
    return(parser)
}

## MAIN

## Parameters
args_p <- get_argparser()
args <- args_p$parse_args(commandArgs(trailingOnly=TRUE))
if(args$verbose){
    write(sprintf('Script summarizing CV prediction results\nGiven arguments:'),stdout())
    write(paste(paste(names(unlist(args)),unlist(args),sep=': '),collapse='\n'),stdout())
}

# source(sprintf("%s/my_utils.R",args$src_path))
source(sprintf("%s/perf_utils.R",args$src_path))

## Merge predictions
pred_mat <- class_num <- NULL
for(pred_file in args$pred_files){
    if(!file.exists(pred_file) | file.info(pred_file)$size==0){
#         stop(sprintf("Empty prediction file %s",pred_file))
        write(sprintf("QUIT: Empty prediction file %s",pred_file), stdout())
        q("no", status=0, runLast=FALSE)
    }
    p <- read.csv(file=pred_file, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    
    if(args$verbose){
        write(sprintf("Prediction file %s: %d rows, %d cols",pred_file,nrow(p),ncol(p)),stdout())
    }
    
    p_fold <- as.numeric( sub( 'fold' , '' , regmatches(pred_file, regexpr('fold[[:digit:]]+',pred_file)) ) )
    if(is.null(pred_mat)){
        pred_mat  <- p
        class_num <- data.frame(fold=p_fold, num0=sum(p$trueClass==0), num1=sum(p$trueClass==1), check.names=FALSE, stringsAsFactors=FALSE)
    }
    else{
        testit::assert(all(colnames(p)==colnames(pred_mat))) # same columns
        testit::assert(length(intersect( pred_mat$sampleID , p$sampleID ))==0) # no sample intersection
        pred_mat  <- rbind(pred_mat,p)
        class_num <- rbind(class_num, c(fold=p_fold, num0=sum(p$trueClass==0), num1=sum(p$trueClass==1)))
    }
}
# Check sample IDs
testit::assert("Sample IDs from all prediction files are not unique.",length(unique(pred_mat$sampleID))==nrow(pred_mat))

## Performance
# Class labels
pred_levels <- sort(unique(pred_mat$trueClass))

# Confusion matrix
pred_cm     <- confusion_mat(true=pred_mat$trueClass, pred=pred_mat$predClass, levels=pred_levels)
# pred_cm_str <- print_confusion_mat(pred_cm)

# Performance
pred_perf   <- all_perf(true=pred_mat$trueClass, pred=pred_mat$predClass, levels=sort(unique(pred_mat$trueClass)), with_names=TRUE)

# AUC
pred_auc    <- c(
    AUC_ROC=get_auc_roc(true=pred_mat$trueClass, pred=1-pred_mat$predCprob, label.ordering=sort(unique(pred_mat$trueClass))),
    AUC_PR =get_auc_pr( true=pred_mat$trueClass, pred=1-pred_mat$predCprob, label.ordering=sort(unique(pred_mat$trueClass)))
)

# Result table
pred_tab  <- c(
    num0=sum(pred_mat$trueClass==0),
    num1=sum(pred_mat$trueClass==1),
    TN=pred_cm$TN, TP=pred_cm$TP, FN=pred_cm$FN, FP=pred_cm$FP,
    pred_perf,
    pred_auc
)
pred_tab  <- matrix(pred_tab, nrow=1, dimnames=list(NULL, names(pred_tab)))

## Plots
if(args$plot){
    pdf(sprintf("%s.pdf",file_path_sans_ext(args$o_file)), width=10, height=5)
    plot_eval(true=pred_mat$trueClass, pred=1-pred_mat$predCprob, label.ordering=sort(unique(pred_mat$trueClass)))
    dev.off()
}

## Save
write.table(x=pred_tab,  file=args$o_file,                                                  sep='\t', row.names=FALSE, quote=FALSE)
write.table(x=class_num, file=sprintf("%s_classnum.csv",file_path_sans_ext(args$o_file)),   sep='\t', row.names=FALSE, quote=FALSE)