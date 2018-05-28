#!/usr/bin/Rscript

plot_pc_var <- function(ev, k=NA){
    ev_df <- data.frame(
        PC=1:length(ev), # principal component
        PCT=100 * ev / sum(ev),
        CumPCT=sapply(1:length(ev), function(i){
            ifelse(i==1, 100*ev[i]/sum(ev), sum(100*ev[1:(i-1)]/sum(ev)) + (100*ev[i]/sum(ev)))
        })
    )
    br <- na.omit(c( k , sapply(0:round(log10(length(ev))), function(x){ 10^x }) ))
    testit::assert(round(max(ev_df$CumPCT))<=100)
    p <- ggplot(ev_df, aes(x=PC, y=CumPCT)) + 
    geom_bar(stat="identity", colour='black', fill='white', aes(width=0.3/(PC))) + 
    scale_x_log10(breaks=br, labels=br) +
    ylab('Cumulative variance proportion [%]') +
    theme(
        panel.grid.major=element_line(colour="#CCCCCC"),
        panel.grid.minor=element_line(colour="white"),
        axis.text=element_text(colour="black"),
        panel.background=element_rect(fill="white"),
        legend.key=element_rect(fill="white")
    )
    if(!is.na(k)){
        p <- p + geom_vline(xintercept=k, linetype="twodash")
    }
    return(p)
}