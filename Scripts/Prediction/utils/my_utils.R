#!/usr/bin/Rscript

## Read in
read_x <- function(x_file, cols_are_factors=FALSE, ordered_factors=FALSE, verbose=TRUE){
    x_mat <- read.csv(file=x_file, sep='\t', header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=TRUE)
    if(cols_are_factors){
        for(i in 1:ncol(x_mat)){
            x_mat[,i] <- factor(x_mat[,i], ordered=ordered_factors)
        }
    }
    if(verbose){
        write(sprintf('x: Read in %s:\n\t> %d rows and %d columns.',x_file,nrow(x_mat),ncol(x_mat)), stdout())
    }
    return(x_mat)
}

read_y <- function(y_file, cols_are_factors=FALSE, ordered_factors=FALSE, verbose=TRUE){
    y_mat <- read.csv(file=y_file, sep='\t', header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=TRUE)
    if(cols_are_factors){
        for(i in 1:ncol(y_mat)){
            y_mat[,i] <- factor(y_mat[,i], ordered=ordered_factors)
        }
    }
    if(verbose){
        write(sprintf('y: Read in %s:\n\t> %d rows and %d columns.', y_file, nrow(y_mat), ncol(y_mat)), stdout())
        write(paste(
            sapply(colnames(y_mat),function(p){
                sprintf("%s: %s",p,paste(
                    sapply(sort(unique(y_mat[,p])),function(v){
                        sprintf("\t%s=%.2f (%d)",as.character(v),100*sum(y_mat[,p]==v,na.rm=TRUE)/nrow(y_mat),sum(y_mat[,p]==v,na.rm=TRUE))
                    }),
                    collapse='; '
                ))
            }),
            collapse='\n'
        ), stdout())
    }
    return(y_mat)
}

read_list <- function(list_file, verbose=TRUE){
    if(is.null(list_file)){ return(NULL) }
    l <- tryCatch(
        as.character(unique(unlist(read.csv(file=list_file, header=FALSE, stringsAsFactors=FALSE, comment.char="#")))),
        error = function(e){ NULL }
    )
    if(!is.null(l)){
        names(l) <- NULL
        if(verbose){
            write(sprintf('List: Read in %s with %d unique entries.',list_file,length(l)), stdout())
        }
    }
    else{
        if(verbose){
            write(sprintf('Could not read list %s',list_file,length(l)), stdout())
        }
    }
    return(l)
}

## Select columns in matrix by names or IDs
sel_mat_cols <- function(mat, incl_cnames=NULL, excl_cnames=NULL, incl_cids=NULL, excl_cids=NULL, verbose=TRUE){
    if(is.null(incl_cnames) & is.null(incl_cids) & is.null(excl_cnames) & is.null(excl_cids)){ return(mat) } # nothing to include/exclude
    testit::assert(
        "Column names AND column IDs given but only one should be given",
        xor(!is.null(incl_cnames)|!is.null(excl_cnames), !is.null(incl_cids)|!is.null(excl_cids))
    )
    if(!is.null(incl_cids) | !is.null(excl_cids)){ # use col names if IDs given
        if(is.null(colnames(mat))){ colnames(mat) <- paste('col',1:ncol(mat),sep='') } # set names if no given
        if(!is.null(incl_cids)){
            testit::assert("Col IDs are outside of the matrix range.", all(incl_cids>=1 & incl_cids<=ncol(mat)))
            incl_cnames <- colnames(mat)[incl_col_ids]
        }
        if(!is.null(excl_cids)){
            testit::assert("col IDs are outside of the matrix range.", all(excl_cids>=1 & excl_cids<=ncol(mat)))
            excl_cnames <- colnames(mat)[excl_col_ids]
        }
    }
    
    keep_cnames <- colnames(mat)
    if(!is.null(incl_cnames)){ keep_cnames <- intersect(keep_cnames,incl_cnames) }
    if(!is.null(excl_cnames)){ keep_cnames <- setdiff(keep_cnames,excl_cnames) }
    testit::assert(all(sapply(keep_cnames, function(x){ x %in% colnames(mat) })))
    write(sprintf(
        '\nMatrix contains %d columns\n\t%d should be included, %d are in matrix\n\t%d should be excluded, %d are in matrix\n\t%d are in include and exclude\n\t%d columns remained',
        ncol(mat),
        ifelse(is.null(incl_cnames),0,length(incl_cnames)), length(intersect(colnames(mat),incl_cnames)),
        ifelse(is.null(excl_cnames),0,length(excl_cnames)), length(intersect(colnames(mat),excl_cnames)),
        length(intersect(excl_cnames,incl_cnames)),
        length(keep_cnames)
    ), stdout())
    return(mat[,keep_cnames,drop=FALSE])
}

## Select rows (samples) in matrix by names or IDs
sel_mat_rows <- function(mat, incl_rnames=NULL, excl_rnames=NULL, incl_rids=NULL, excl_rids=NULL, verbose=TRUE){
    if(is.null(incl_rnames) & is.null(incl_rids) & is.null(excl_rnames) & is.null(excl_rids)){ return(mat) } # nothing to include/exclude
    testit::assert(
        "Row names AND row IDs given but only one should be given",
        xor(!is.null(incl_rnames)|!is.null(excl_rnames), !is.null(incl_rids)|!is.null(excl_rids))
    )
    if(!is.null(incl_rids) | !is.null(excl_rids)){ # use row names if IDs given
        if(is.null(rownames(mat))){ rownames(mat) <- paste('row',1:nrow(mat),sep='') } # set names if no given
        if(!is.null(incl_rids)){
            testit::assert("Row IDs are outside of the matrix range.", all(incl_rids>=1 & incl_rids<=nrow(mat)))
            incl_rnames <- rownames(mat)[incl_col_ids]
        }
        if(!is.null(excl_rids)){
            testit::assert("Row IDs are outside of the matrix range.", all(excl_rids>=1 & excl_rids<=nrow(mat)))
            excl_rnames <- rownames(mat)[excl_col_ids]
        }
    }
    
    keep_rnames <- rownames(mat)
    if(!is.null(incl_rnames)){ keep_rnames <- intersect(keep_rnames,incl_rnames) }
    if(!is.null(excl_rnames)){ keep_rnames <- setdiff(keep_rnames,excl_rnames) }
    testit::assert(all(sapply(keep_rnames, function(x){ x %in% rownames(mat) })))
    write(sprintf(
        '\nMatrix contains %d rows (samples)\n\t%d should be included, %d are in matrix\n\t%d should be excluded, %d are in matrix\n\t%d are in include and exclude\n\t%d rows remained',
        nrow(mat),
        ifelse(is.null(incl_rnames),0,length(incl_rnames)), length(intersect(rownames(mat),incl_rnames)),
        ifelse(is.null(excl_rnames),0,length(excl_rnames)), length(intersect(rownames(mat),excl_rnames)),
        length(intersect(excl_rnames,incl_rnames)),
        length(keep_rnames)
    ), stdout())
    return(mat[keep_rnames,,drop=FALSE])
}

## Consistency w.r.t. rownames/samples
make_same_rows <- function(mats, verbose=TRUE){
    if(verbose){ write("Making matrices consistent w.r.t. row names: same names, same order.", stdout()) }
    common <- rownames(mats[[1]])
    for(i in 1:length(mats)){
        testit::assert(sprintf("\tMatrix %s has no row names", ifelse(is.null(names(mats)),as.character(i),names(mats)[i])), !is.null(rownames(mats[[i]])))
        common <- intersect(common, rownames(mats[[i]]))
        if(verbose){
            write(sprintf("\tMatrix %s has %d rows", ifelse(is.null(names(mats)),as.character(i),names(mats)[i]), nrow(mats[[i]])), stdout())
        }
    }
    common <- sort(common)
    if(verbose){ write(sprintf("There are %d row names contained in all matrices.", length(common)), stdout()) }
    for(i in 1:length(mats)){
        testit::assert(all(sapply(common, function(x){ x %in% rownames(mats[[i]]) })))
        mats[[i]] <- (mats[[i]])[common,]
    }
    return(mats)
}

## Filter matrix by columns: empty and constant ones are removed in any case
filter_mat_cols <- function(mat, missing_value=NA, max_fq=100, max_missfq=0, cores=2, verbose=TRUE){
    # if missing are not NA, set to NA
    if(!is.na(missing_value)){
        mat_missing <- (mat==missing_value)
        testit::assert(sprintf("Missing value is \"%s\" but there are NAs in the matrix.",missing_value),all(!is.na(mat)))
        mat[mat==missing_value] <- NA
        testit::assert(all(mat_missing==is.na(mat)))
        if(verbose){ write(sprintf('IMPORTANT: Missing value \"%s\" was replaced ba NA',missing_value),stdout()) }
    }
    
    library(parallel)
    successed <- FALSE; tried <- 0
    while(!successed & tried <5){
        tryCatch({
            cl <- makeCluster(cores, type = "FORK")
#             clusterExport(cl=cl, varlist=c('mat','max_fq','max_missfq'), envir=environment())
            keep <- parSapply(cl=cl, X=colnames(mat), FUN=function(f){ 
                x <- as.vector(mat[,f])
                if(all(is.na(x))){ return(FALSE) } # all missing -> do not keep
                x_tab <- table(x[!is.na(x)]) # table without missing
                x_mm  <- length(unique(x[!is.na(x)]))==1 # TRUE if monomorphic (only one non-missing value)
                x_mf  <- 100*max(x_tab)/sum(x_tab) # maximal genotype freq. [%]
                x_mr  <- 100*sum(is.na(x))/length(x) # missing rate [%]
                # keep if: non-monomorphic AND max. freq. <= t AND missing rate <= t
                return( (!x_mm) & x_mf<=max_fq & x_mr<=max_missfq )
            }, USE.NAMES=FALSE)
            testit::assert(is.vector(keep))
            testit::assert(length(keep) == ncol(mat))
            stopCluster(cl)
            successed <- TRUE
        }, error=function(e){ tried <- tried + 1 }
        )
    }
    if(!successed){ stop(sprintf("Could not perform parallel computing")) }
    
    if(verbose){
        write(sprintf(
            "Matrix column filtering: From %d columns, %d were removed (monomorphic, max. freq. rate > %0.3f, missing rate > %0.3f): %d remained",
            ncol(mat), ncol(mat)-sum(keep), max_fq, max_missfq, sum(keep)
        ),stdout())
    }
    return(mat[,keep,drop=FALSE])
}

## Remove correlated columns by hierarchical clustering
# Get medoid = cluster member with smallest mean distance within the cluster
get_medoid <- function(d_dist_mat, clusters, cluster){
    # names of cluster members
    members <- names(clusters)[clusters==cluster]
    if(length(members)==1){return(members)}
    # mean distance per member within the cluster -> vector of mean dist. for each member
    ave_dist <- apply(d_dist_mat[members,members], 1, mean)
    # medoid = member with minimal mean distance -> name of member with min. mean dist.
    medoid <- names(ave_dist)[which.min(ave_dist)]
    return(medoid)
}
# Apply correlation based hierarchical clustering to remove correlated columns/predictors
remove_corr_hclust <- function(mat, mcor, cut_h=1-(0.7**2), cut_k=10, use_h=TRUE, aggl.method='average', verbose=TRUE){
    # distance
    mcor_dist_mat <- 1-(mcor**2) ## dist = 1 - r^2 ref: http://wpicr.wpic.pitt.edu/WPICCompGen/hclust/hclust.htm
    testit::assert(all(!is.na(mcor_dist_mat)))
    mcor_dist <- as.dist(mcor_dist_mat, diag=TRUE, upper=TRUE)
    # hclust
    cor_tree <- hclust(d=mcor_dist, method=aggl.method)
    cor_tree$height <- round(cor_tree$height,6) # to avoid "the 'height' component of 'tree' is not sorted (increasingly)" (https://stat.ethz.ch/pipermail/r-help/2008-May/163409.html)
    # cut the tree -> vectors of cluter IDs with names of clustered elements
    cor_clust <- NULL
    if(use_h){
        cor_clust <- cutree(tree=cor_tree, h=cut_h)
    }
    else{
        cor_clust <- cutree(tree=cor_tree, k=cut_k)
    }
    # get cluster representatives (medoids) for each cluster
    clust_medoids <- sapply(sort(unique(cor_clust)), function(x){
        get_medoid(d_dist_mat=mcor_dist_mat, clusters=cor_clust, cluster=x)
    })
#     clust_medoids <- unlist(clust_medoids)
    # which should be removed
    to_remove <- rep(TRUE, ncol(mat)); names(to_remove) <- colnames(mat)
    to_remove[clust_medoids] <- FALSE # keep medoids
    if(verbose){
        write(sprintf(
            'Removed corr. columns using hclust (dist=1-cor^2, aggl. method %s, cut at %s = %.5f): From %d col.s %d were kept (%d removed).',
             aggl.method, ifelse(use_h,'h','k'), ifelse(use_h,cut_h,cut_k),
             ncol(mat),sum(!to_remove),sum(to_remove)
        ), stdout())
    }
    # return matrix without correlated columns and column clusters
    return(list(
        mat=mat[,!to_remove,drop=FALSE],
        mat_cl=data.frame(ID=names(cor_clust),Cluster=cor_clust,stringsAsFactors=FALSE)
    ))
}

## Printing
# for printing with sprintf
# vec_as_string <- function(vec, prefix='\t', print_width=5){
#     if(length(vec)==0){ return('') }
#     # split into chuncks
#     vec_split <- split(vec, ceiling(seq_along(1:length(vec))/print_width))
#     # collapse each chunck
#     vec_str <- lapply(vec_split,function(x){ x[1]<-paste(prefix,x[1],sep=''); return( paste(x,collapse=', ') ) })
#     # collapse chuncks together
#     vec_str <- paste(unlist(vec_str),collapse=',\n')
#     return(vec_str)
# }