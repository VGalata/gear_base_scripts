#!/usr/bin/Rscript

## F = false predictions
## T = true predictions
## TP = true resistant
## TN = true not resistant
## Error rate: err = (FP+FN) / (TP+TN+FP+FN) = F / all
## Accuracy: acc = (TP+TN) / (TP+TN+FP+FN) = T / all
## Sensitivity: sens = Recall = TPR = TP / (TP+FN)
## Specificity: spec = TNR = TN / (TN+FP)
## Precision: prec = TP / (TP+FP)
## FPR = FP / (FP+TN)
## FNR = FN / (FN+TP)
## Fmeasure: Fmeasure = 2 * (Recall * Precision) / (Recall + Precision)
## Geometric mean with Recall Specificity: gm_RS = sqrt( Recall * Specificity )
## Geometric mean with Recall Precision: gm_RP = sqrt( Recall * Precision )

## PACKAGES
suppressMessages(require(ROCR)) # computing prediction and performance objects, plotting
suppressMessages(require(MESS)) # compute AUC for a given curve



#############################
## CLASSIFICATION, CLASSES ##
#############################
## Note: prediction contains class labels

## Compute confusion matrix and TN, TP, FN, FP
## Note: Factors in true and predicted class labels: Positive label > Negative label, e.g. (TRUE/FALSE, 1/-1, 1/0)
confusion_mat <- function(true, pred, levels){
    results <- list(conf=NULL, TN=0, TP=0, FN=0, FP=0)
    
    results$conf <- table(true=factor(true,levels=levels), pred=factor(pred,levels=levels))
    
    results$TP <- results$conf[2,2]
    results$TN <- results$conf[1,1]
    results$FP <- results$conf[1,2]
    results$FN <- results$conf[2,1]
    
    return(results)
}

print_confusion_mat <- function(conf_m){
    conf_m <- conf_m$conf # get table object
    # paste: <row v.>,<col v.>:<cell count> for each cell sep. by "|"
    conf_m_str <- paste(
        paste(
            apply(expand.grid(rownames(conf_m),colnames(conf_m)),1,function(s){ paste(s,collapse=',') }),
            as.vector(conf_m),
            sep=':'
        ),
        collapse=' | '
    )
    # add names: here true vs. pred
    conf_m_str <- paste(paste(names(dimnames(conf_m)),collapse=' vs '),conf_m_str,sep=' : ')
    return(conf_m_str)
}

## Error rate = (FP+FN) / (TP+TN+FP+FN) = F / all
err <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round((conf_m$FP+conf_m$FN)/sum(conf_m$conf),3))
}

## Accuracy = (TP+TN) / (TP+TN+FP+FN) = T / all
acc <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round((conf_m$TP+conf_m$TN)/sum(conf_m$conf),3))
}

## Balanced Accuracy = (Sens+Spec) / 2
b_acc <- function(true, pred, levels){
    return((sens(true, pred, levels)+spec(true, pred, levels))/2)
}

## Sensitivity = Recall = TPR = TP / (TP+FN)
sens <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round(conf_m$TP/(conf_m$TP+conf_m$FN),3))
}

## Specificity = TNR = TN / (TN+FP)
spec <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round(conf_m$TN/(conf_m$TN+conf_m$FP),3))
}

## Precision = TP / (TP+FP)
prec <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round(conf_m$TP/(conf_m$TP+conf_m$FP),3))
}

## NPV = TN / (TN + FN)
NPV <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round(conf_m$TN/(conf_m$TN+conf_m$FN),3))
}

## FPR = FP / (FP+TN)
fpr <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round(conf_m$FP/(conf_m$FP+conf_m$TN),3))
}

## FNR = FN / (FN+TP)
fnr <- function(true, pred, levels){
    conf_m <- confusion_mat(true=true, pred=pred, levels=levels)
    return(round(conf_m$FN/(conf_m$FN+conf_m$TP),3))
}

## F-measure = 2 * (Recall * Precision) / (Recall + Precision)
Fmeasure <- function(true, pred, levels){
    p <- prec(true=true, pred=pred, levels=levels)
    r <- sens(true=true, pred=pred, levels=levels)
    return(round(2*(r*p)/(r+p),3))
}

## Geometric mean RS = sqrt( Recall * Specificity )
gm_RS <- function(true, pred, levels){
    s <- spec(true=true, pred=pred, levels=levels)
    r <- sens(true=true, pred=pred, levels=levels)
    return(round(sqrt(s*r),3))
}

## Geometric mean RP = sqrt( Recall * Precision )
gm_RP <- function(true, pred, levels){
    p <- prec(true=true, pred=pred, levels=levels)
    r <- sens(true=true, pred=pred, levels=levels)
    return(round(sqrt(r*p),3))
}

## all above measures
all_perf <- function(true, pred, levels, with_names=FALSE){
    res <-c(    err(true=true, pred=pred, levels=levels),
                acc(true=true, pred=pred, levels=levels),
                b_acc(true=true, pred=pred, levels=levels),
                sens(true=true, pred=pred, levels=levels),
                spec(true=true, pred=pred, levels=levels),
                prec(true=true, pred=pred, levels=levels),
                NPV(true=true, pred=pred, levels=levels),
                fpr(true=true, pred=pred, levels=levels),
                fnr(true=true, pred=pred, levels=levels),
                Fmeasure(true=true, pred=pred, levels=levels),
                gm_RS(true=true, pred=pred, levels=levels),
                gm_RP(true=true, pred=pred, levels=levels)
            )
#     if(with_names){ names(res) <- c('err','acc','b_acc','sens','spec','prec','NPV','fpr','fnr','Fmeasure','gm_RS','gm_RP')}
    if(with_names){ names(res) <- c('ERR','ACC','B_ACC','SENS','SPEC','PREC','NPV','FPR','FNR','Fmeasure','gm_RS','gm_RP')}
    return(res)
}



###########################
## CLASSIFICATION, PROBS ##
###########################
## Note: prediction contains class probabilities

## Create the prediction object from true graph G and edge conf. level matrix eGprob and the prior
get_pred <- function(true, pred, label.ordering=NULL){
    return(prediction(predictions=pred, labels=true, label.ordering=label.ordering))
}

## Given performance object, compute the AUC (NaN, infinity are set to 1)
compute_auc <- function(perf){
    x <- perf@x.values[[1]]; y <- perf@y.values[[1]]
    # NaN and infinity => 1 (bit return FALSE for is.finite(.))
    x[!is.finite(x)] <- 1; y[!is.finite(y)] <- 1
    # return the AUC
    return(abs(round(MESS::auc(x=x,y=y),3)))
}

## Given the prediction object, get AUC of the precision-recall curve
get_auc_pr <- function(true, pred, label.ordering=NULL){
    pred_o <- get_pred(true=true, pred=pred, label.ordering=label.ordering)
    return(compute_auc(performance(pred_o, measure="prec", x.measure="rec")))
}

## Given the prediction object, get AUC of the TPR-FPR curve
get_auc_roc <- function(true, pred, label.ordering=NULL){
    pred_o <- get_pred(true=true, pred=pred, label.ordering=label.ordering)
    return(abs(round(performance(pred_o, measure="auc")@y.values[[1]],3)))
}


## Given pred. object, draw the precision-recall curve
plot_prec_recall <- function(true, pred, info='', label.ordering=NULL){
    pred_o <- get_pred(true=true, pred=pred, label.ordering=label.ordering)
    if(is.null(pred_o)){return()}
    perf <- performance(pred_o, measure="prec", x.measure="rec")
    auc.value <- compute_auc(perf)
    plot(perf,main=paste(info,'AUC=',auc.value,sep=''), xlim=c(0,1), ylim=c(0,1)) # plot the curve
    # AUC
    x <- perf@x.values[[1]]; y <- perf@y.values[[1]]; x[!is.finite(x)] <- 1; y[!is.finite(y)] <- 1
    polygon(x=c(0,x,1),y=c(0,y,0),col='gray')
}

## Given pred. object, draw the TPR-FPR curve
plot_roc <- function(true, pred, info='', label.ordering=NULL){
    pred_o <- get_pred(true=true, pred=pred, label.ordering=label.ordering)
    if(is.null(pred_o)){return()}
    perf <- performance(pred_o, measure="tpr", x.measure="fpr")
    auc.value <- abs(round(performance(pred_o, measure="auc")@y.values[[1]],3))
    plot(perf,main=paste(info,'AUC=',auc.value,sep=''), xlim=c(0,1), ylim=c(0,1)) # plot the curve
    # AUC
    x <- perf@x.values[[1]]; y <- perf@y.values[[1]]; x[!is.finite(x)] <- 1; y[!is.finite(y)] <- 1
    polygon(x=c(x,1),y=c(y,0),col='gray')
    # base line
    lines(x=seq(0,1,0.1),y=seq(0,1,0.1),col='blue')
}

## Given pred. object, draw the TPR-FPR curve and the PR-curve in one plot
plot_eval <- function(true, pred, info='', label.ordering=NULL){
    par(mfrow=c(1,2))
    plot_roc(true=true, pred=pred, info=info, label.ordering=label.ordering)
    plot_prec_recall(true=true, pred=pred, info=info, label.ordering=label.ordering)
}

## NEW BEGIN

## Get ROC curve coordinates
get_roc <- function(true, pred, info='', label.ordering=NULL){
    pred_o <- get_pred(true=true, pred=pred, label.ordering=label.ordering)
    if(is.null(pred_o)){return(NULL)}
    perf <- performance(pred_o, measure="tpr", x.measure="fpr")
    x <- perf@x.values[[1]]; y <- perf@y.values[[1]]; x[!is.finite(x)] <- 1; y[!is.finite(y)] <- 1
#     print(get_auc_roc(true,pred,label.ordering))
    return(data.frame(x=x,y=y))
}
## Draw ROC curves for multiple results (runs) and the ROC curve of the average results
plot_ave_roc <- function(y, pred_probs, main){
    # compute ROC curves for each CV run
    rocs <- lapply(1:ncol(pred_probs), function(s){return(get_roc(true=y, pred=pred_probs[,s], label.ordering=levels(y)))})
    # plot the average ROC curve
    plot_roc(true=y, pred=apply(pred_probs,1,mean), info=main, label.ordering=levels(y))
    # plot all the other ROC curves
    test<-lapply(seq_along(rocs),function(i) lines(rocs[[i]]$x,rocs[[i]]$y,col='black',lty=3))
    # add legend
    legend(x="bottomright", lty=c(1,3,1), col=c('black','black','blue'), legend=c('mean class prob.s','CV repetitions','baseline'), inset=0.05)
}

## NEW END


################
## REGRESSION ##
################
## Residual sum of squares
RSS <- function(y,y_pred){
    sum((y-y_pred)**2)
}

## Mean squared error
MSE <- function(y,y_pred){
    RSS(y,y_pred)/length(y)
}

## Correlation coefficient
corr_coeff <- function(y, y_pred){
    cor(y,y_pred, method='pearson')
}


################
## 2x2 TABLES ##
################
# Check whether top-left/bottom-right diagonal has more samples then bottom-left/top right
diag_check <- function(tab){
    return((tab[1,1]+tab[2,2])>=(tab[1,2]+tab[2,1]))
}
table2x2_perf <- function(tables, x_mat, y_mat, verbose=TRUE){
    if(verbose){
        print('Computing single feature classification performance...')
        print('If 1st x-value->1st y-value and 2nd x-value->2nd y-value then using:')
        print(matrix(c('TN','FN','FP','TP'),nrow=2,ncol=2))
        print('else:')
        print(matrix(c('FN','TN','TP','FP'),nrow=2,ncol=2))
    }
    mat <- data.frame(  acc= unlist(lapply(tables, function(cm){ifelse(diag_check(cm),sum(diag(cm))/sum(cm),(cm[1,2]+cm[2,1])/sum(cm))})),
                        sens=unlist(lapply(tables, function(cm){ifelse(diag_check(cm),cm[2,2]/sum(cm[2,]),cm[1,2]/sum(cm[1,]))})),
                        spec=unlist(lapply(tables, function(cm){ifelse(diag_check(cm),cm[1,1]/sum(cm[1,]),cm[2,1]/sum(cm[2,]))})),
                        prec=unlist(lapply(tables, function(cm){ifelse(diag_check(cm),cm[2,2]/sum(cm[,2]),cm[1,2]/sum(cm[,2]))})),
                        stringsAsFactors=FALSE, check.names=FALSE, row.names=names(tables)
                        )
    return(mat)
}