#!/usr/bin/Rscript

write("Rscript: resistance-date tests\n", stdout())

## LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(reshape2))
suppressMessages(library(exactRankTests))

## TEST
# args <- list(
#     'data'='~/git_repos/Bacteria/Data/10k_bacs_class_res/data/data.tsv',
#     'obanme'='',
#     'min_num'=50,
#     # 'test_min_total'=8,
#     'test_min_group'=10,
#     'adjust'='fdr',
#     'alpha'=0.01,
#     'alt'='two.sided',
#     'src_path'='/home/vgalata/git_repos/Bacteria/SubProject_ClassStat/Pipeline/scripts'
# )

## Command line argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--data', required=TRUE)
    parser$add_argument('--ords', required=TRUE)
    parser$add_argument('--min_num', default=50, type="integer")
    # parser$add_argument('--test_min_total', default=4, type="integer")
    parser$add_argument('--test_min_group', default=10, type="integer")
    parser$add_argument('--adjust', choices=p.adjust.methods)
    parser$add_argument('--alpha', default=0.05, type="double")
    parser$add_argument('--alt', choices=c('two.sided', 'less', 'greater'), required=TRUE)
    parser$add_argument('--src_path', help='')
    return(parser)
}

## FILES & DATA
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## Help functions
suppressMessages(source(sprintf('%s/utils/utils.R',args$src_path)))

## Log file
log_file <- sprintf('%s.log', args$ords)
write("resistance-date tests\n", file=log_file, append=FALSE)
save_session_info(log_file)
save_args(args, log_file)

## Data
tabs <- split_data(tab=read_data(file_name=args$data))

## Samples
# filter
ids <- sapply(rownames(tabs$mics), function(s_id){
    !is.na(tabs$dates[s_id, 'Date']) & # has date
    !all(is.na(tabs$res[s_id,])) # has MIC profile
})
# get IDs
ids <- rownames(tabs$mics)[ids]
# log
write(sprintf("Filtering by date/resistance profile: %d samples", length(ids)), file=log_file, append=TRUE)

## Test data
test_data <- data.frame(
    Date=tabs$dates[ids,'Date'],
    Year=sapply(tabs$dates[ids,'Date'], function(x){ as.numeric(format(x, '%Y')) }), # extract year from date
    species=convert_taxa(tabs$tax[ids,'species']),
    tabs$res[ids,],
    row.names=ids, check.names=FALSE, stringsAsFactors=FALSE
)
# discard taxa
test_data <- discard_taxa(df=test_data, min_num=args$min_num)
# remove samples with taxon == "Rest"
test_data <- test_data[test_data$species != "Rest",]
# unique (remaining) taxa
taxa <- sort(unique(test_data$species))
# log
write(sprintf("Filtering by number of samples per taxon: %d samples, %d taxa", nrow(test_data), length(taxa)), file=log_file, append=TRUE)

## Test
df_rows <- taxa
df_cols <- sort(setdiff(colnames(test_data), c('Date', 'Year', 'species')))

tmp <- data.frame(matrix(NA, nrow=length(df_rows), ncol=length(df_cols), dimnames=list(df_rows, df_cols)), check.names=FALSE)

test_results <- list(
    test_data=test_data,
    skipped=tmp, # test skipped because ...
    pvalues=tmp, # test p-values
    pvalues_adj=tmp, # adjusted p-values
    tests=list() # test object
)

# run test
pb <- txtProgressBar(min=0, max=length(df_rows), initial=0, style=3)
for(taxon in df_rows){
    # subset
    df <- subset(test_data, species==taxon)

    # checks
    if(nrow(df)<=2){ # only two samples
        test_results$skipped[taxon,] <- rep('TWO_SAMPLES', length(df_cols))
        next
    }
    if(length(unique(df$Year))==1){ # only one unique year value
        test_results$skipped[taxon,] <- rep('ONE_DATE', length(df_cols))
        next
    }

    # test for each drug
    for(drug in df_cols){
        # split samples into two groups
        groupA <- df[,drug] != 'R'
        groupB <- df[,drug] == 'R'

        # checks
        if(all(is.na(df[,drug]))){
            test_results$skipped[taxon,drug] <- 'NA_ALL'
            next
        }
        testit::assert(!any(is.na(df[,drug])))
        # if(nrow(df) < args$test_min_total | sum(groupA) < args$test_min_group | sum(groupB) < args$test_min_group){
        if(sum(groupA) < args$test_min_group | sum(groupB) < args$test_min_group){
            # if(nrow(df) < args$test_min_total){
            #     test_results$skipped[taxon,drug] <- 'NUM_ALL'
            # } else
            if(sum(groupA) < args$test_min_group){
                test_results$skipped[taxon,drug] <- 'NUM_A'
            } else if(sum(groupB) < args$test_min_group){
                test_results$skipped[taxon,drug] <- 'NUM_B'
            }
            next
        }
        if(length(unique(df[,'Year']))==1){
            test_results$skipped[taxon,drug] <- 'ALL_SAME_YEAR'
            next
        }
        if(length(unique(df[groupA,'Year']))==1 & length(unique(df[groupA,'Year']))==1){
            test_results$skipped[taxon,drug] <- 'ALL_SAME_YEAR_PER_GROUP'
            next
        }
        if(length(unique(df[groupA,'Year']))==1){
            test_results$skipped[taxon,drug] <- 'ALL_SAME_YEAR_GROUP_A'
            next
        }
        if(length(unique(df[groupB,'Year']))==1){
            test_results$skipped[taxon,drug] <- 'ALL_SAME_YEAR_GROUP_B'
            next
        }

        # test
        set.seed(42) # used exactRankTests::wilcox.exact call is not deterministic
        test_test <- exactRankTests::wilcox.exact(
            x=df[groupA,'Year'],
            y=df[groupB,'Year'],
            paired=FALSE,
            alternative=args$alt,
            conf.int=TRUE, conf.level=1-args$alpha
        )
        testit::assert(! 'MY_TOTAL' %in% names(test_test))
        testit::assert(! 'MY_GROUP_A' %in% names(test_test))
        testit::assert(! 'MY_GROUP_B' %in% names(test_test))
        test_test$MY_TOTAL <- nrow(df)
        test_test$MY_GROUP_A <- groupA
        test_test$MY_GROUP_B <- groupB
        test_test$MY_GROUP_A_SIZE <- sum(groupA)
        test_test$MY_GROUP_B_SIZE <- sum(groupB)

        # save
        test_results$tests <- c(test_results$tests, list(test_test))
        names(test_results$tests)[length(test_results$tests)] <- sprintf("%s:%s", taxon, drug)
        test_results$pvalues[taxon, drug] <- test_test$p.value
    }
    setTxtProgressBar(pb, which(df_rows==taxon))
}
close(pb)

## Adjust p-values
# reshape p-value matrix: taxon, drug, p-value
pv <- melt(cbind(Taxon=rownames(test_results$pvalues), test_results$pvalues), id.vars="Taxon")
# adjust
pv$adj <- p.adjust(pv$value, method=args$adjust)
# save in taxon x drug matrix
for(taxon in df_rows){
    for(drug in df_cols){
        test_results$pvalues_adj[taxon, drug] <- pv[pv$Taxon==taxon & pv$variable==drug,'adj']
    }
}

## Save
ofile <- args$ords
write(sprintf("Saving to %s", ofile), file=log_file, append=TRUE)
saveRDS(object=test_results, file=ofile, compress=TRUE)
