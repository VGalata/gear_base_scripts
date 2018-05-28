#!/usr/bin/Rscript

# http://www.gettinggeneticsdone.com/2010/07/qq-plots-of-p-values-in-r-using-base.html
# http://GettingGeneticsDone.blogspot.com/
# See http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Define the function
ggd.qqplot = function(pvector, main=NULL, ...) {
    pvector <- as.vector(na.omit(pvector)) ## V. G. 2016.02.19
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o)))
    lines(e,e,col="red")
}