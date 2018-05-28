#!/usr/bin/Rscript

## LIBS
suppressMessages(library(argparse))
# suppressMessages(library(testit))
# suppressMessages(library(tools)) # file name (without ext.) from path etc.
suppressMessages(library(rpart)) # trees
# suppressMessages(library(rpart.plot)) # for plotting

## FUNCTIONS
# Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    # required parameters
    parser$add_argument('-x_file',  help='Feature file', required=TRUE)
    parser$add_argument('-mod_file',  help='Feature file', required=TRUE)
    parser$add_argument('-o_file',   help='Output directory', required=TRUE)
    
    parser$add_argument('-y_file',  help='')
    parser$add_argument('-y_name',  help='')
    parser$add_argument('-samples',  help='')
    parser$add_argument('--transpose_x', action='store_true')
    
    parser$add_argument('--src_path',   help='Print additional information.', default='.')
    
    parser$add_argument('--verbose',    help='Print additional information.', action='store_true')
    
    return(parser)
}

## MAIN
## Args
args    <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

if(args$verbose){
    write(sprintf('\n\n'), stdout())
}

## Model
pred_mod <- readRDS(file=args$mod_file)

if( length(setdiff( sapply(labels(pred_mod), function(x){ strsplit(x,'=')[[1]][1] }), c('root'))) == 0 ){
    if(args$verbose){
        write(sprintf('empty model %s',args$mod_file), stdout())
    }
    q("no", status=0, runLast=FALSE)
}

## Data
# X
# Read
x       <- read.csv(file=args$x_file, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
if(args$verbose){
    write(sprintf('Read in x: %d rows and %d columns',nrow(x),ncol(x)), stdout())
}
# Transpose
if(args$transpose_x){
    x <- data.frame(t(x), check.names=FALSE)
}
if(args$verbose){
    write(sprintf('Transposed x: %d rows and %d columns',nrow(x),ncol(x)), stdout())
}
# Select features contained in model (to reduce dimensions)
x <- x[,names(attr(pred_mod,'xlevels')),drop=FALSE]
# Factors
for(i in 1:ncol(x)){ x[,i] <- factor(x[,i], levels=c('0','1')) }
if(args$verbose){
    write(sprintf('Casted col.s in x to factors'), stdout())
}

# Y
y <- NULL
if(!is.null(args$y_file)){
    # Read
    y <- read.csv(file=args$y_file, header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
    if(args$verbose){
        write(sprintf('Read in y: %d rows and %d columns',nrow(y),ncol(y)), stdout())
    }
    # Factors
    for(i in 1:ncol(y)){ y[,i] <- factor(y[,i], levels=c('0','1')) }
    # Select column by name
    if(args$verbose){
        write(sprintf('Read in x: %d rows and %d columns',nrow(x),ncol(x)), stdout())
    }
    testit::assert(!is.null(args$y_name))
    args$y_name <- gsub('_','/',args$y_name) # replace "_" by "/" in name
    testit::assert(sprintf("y column name %s not in y",args$y_name),args$y_name %in% colnames(y))
    y <- y[,args$y_name,drop=FALSE]
}

# Samples
samples <- NULL
if(!is.null(args$samples)){
    samples <- as.character( unlist( read.csv(file=args$samples, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE) ) )
    if(args$verbose){
        write(sprintf('Read in samples: %d',length(samples)), stdout())
    }
    # select those contained in x and y with non-NA y
    samples <- samples[ sapply(samples, function(s){ s %in% rownames(x) & s %in% rownames(y) & !is.na(unlist(y[s,])) }) ]
    if(args$verbose){
        write(sprintf('Selecting y column %s',args$y_name), stdout())
    }
    x <- x[samples,,drop=FALSE]
    y <- y[samples,,drop=FALSE]
}

##
# http://stackoverflow.com/questions/22315394/factor-has-new-levels-error-for-variable-im-not-using
# for (var_name in names(attr(pred_mod,'xlevels'))){
#     if(! all(attr(pred_mod,'xlevels')[[var_name]] == c('0','1'))){
#         attr(pred_mod,'xlevels')[[var_name]] <- c('0','1')
#     }
# }
# print(attr(pred_mod,'xlevels')[['var_1805']])

## Predict
pred_mat <- data.frame(
    sampleID=rownames(x),
    predClass= predict(pred_mod, x, type='class'),
    predCprob=(predict(pred_mod, x, type='prob'))[,'0'],
    check.names=FALSE
)
write.table(x=pred_mat, file=args$o_file, row.names=FALSE, quote=FALSE, sep='\t')

## Perf (if possible)
if(!is.null(y)){
    source(sprintf("%s/utils/perf_utils.R",args$src_path))
    # Class labels
    pred_levels <- c('0','1')
    # Add true class
    pred_mat$trueClass <- unlist(y)
    # Confusion matrix
    pred_cm     <- confusion_mat(true=unlist(y), pred=pred_mat$predClass, levels=pred_levels)
    if(args$verbose){
        write(sprintf('Confusion matrix:'), stdout())
        print(pred_cm$conf)
    }
    # Performance
    pred_perf   <- all_perf(true=pred_mat$trueClass, pred=pred_mat$predClass, levels=pred_levels, with_names=TRUE)
    # AUC
    pred_auc    <- c(
        AUC_ROC=get_auc_roc(true=pred_mat$trueClass, pred=1-pred_mat$predCprob, label.ordering=pred_levels),
        AUC_PR =get_auc_pr( true=pred_mat$trueClass, pred=1-pred_mat$predCprob, label.ordering=pred_levels)
    )

    # Result table
    pred_tab  <- c(
        pheno=args$y_name,
        num0=sum(pred_mat$trueClass==0),
        num1=sum(pred_mat$trueClass==1),
        TN=pred_cm$TN, TP=pred_cm$TP, FN=pred_cm$FN, FP=pred_cm$FP,
        pred_perf,
        pred_auc
    )
    pred_tab  <- matrix(pred_tab, nrow=1, dimnames=list(NULL, names(pred_tab)))
    if(args$verbose){
        write(sprintf('Preformance:'), stdout())
        print(pred_tab)
    }
    write.table(x=pred_tab, file=sprintf('%s.perf',args$o_file), row.names=FALSE, quote=FALSE, sep='\t')
}
