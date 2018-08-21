#' Conditional Local False Discovery Rate
#'
#' @description
#' An implementaiton of the conditional local false discovery rate procedure of
#' Cai and Sun (2009). This is a near copy of the implementation provided in the
#' \pkg{IHWpaper} package, \code{\link[IHWpaper]{clfdr}} at
#' \href{https://github.com/nignatiadis/IHWpaper/blob/c2f6799/R/stratified_methods.R}{the corresponding GitHub repo},
#' written by Nikos Ignatiadis and made available under the GPL-3 license.
#'
#' The original code has been modified to only return a vector of adjusted
#' p-values and check the input of the \code{lfdr_estimation=} parameter.  
#' 
#' @param p numeric vector of unadjusted p-values.
#' @param groups factor to which different hypotheses belong.
#' @param lfdr_estimation method used to estimate the local fdr; must be
#'        one of "fdrtool" or "locfdr". (default = "fdrtool")
#'
#' @details
#' If the covariate or variable for stratifying p-values is continuous, it should
#' be converted into a discrete set of "groups" before applying this function.
#' This can be accomplished by passing the continuous variable to \code{IHW::groups_by_filter},
#' \code{ggplot2::cut_number}, or any similar function. See examples for
#' uses of these functions.
#' 
#' @return
#' Numeric vector of adjusted p-values of equal length and order
#' as the to input vector of p-values.
#'
#' @references
#' Reference for the conditional local false discovery rate (underlying theory):
#' 
#' Cai, T.T., Sun, W. (2009) Simultaneous testing of grouped hypotheses: Finding needles in multiple haystacks. Journal of the American Statistical Association, 104(488):1467-1481. \url{https://doi.org/10.1198/jasa.2009.tm08415}
#' 
#' Reference for the original implementation (underlying code):
#' 
#' Ignatiadis, N. (2017) IHWpaper: Reproduce figures in IHW paper. R package version 1.7.0. \url{https://doi.org/doi:10.18129/B9.bioc.IHWpaper}
#' 
#' @examples
#' ## generate example p-values and covariate values
#' pv <- runif(100)
#' x <- sample(100)
#'
#' ## discretize 'x' into 10 bins
#' xd <- IHW::groups_by_filter(xd, nbins = 10)
#'
#' ## adjust p-values, grouping by covariate bin
#' adjpv <- clfdr(p = pv, groups = xd)
#' 
#' @importFrom fdrtool fdrtool
#' @importFrom locfdr locfdr
#' @export
#' @author Patrick Kimes
clfdr <- function(p, groups, lfdr_estimation = c("fdrtool", "locfdr")) {

    lfdr_estimation <- match.arg(lfdr_estimation)
    
    ## estimate local fdr within each stratum first
    lfdr_res <- lfdr_fit(p, groups, lfdr_estimation = lfdr_estimation)
    lfdrs <- lfdr_res$lfdr
    
    ## Remark:
    ## When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
    ## we get monotonic adjusted p-values as a function of the p-values
    ## This is mainly needed for grenander based lfdrs, with most other
    ## lfdr estimation methods lfdr ties are not a problem usually
    
    ## now use the rejection rule described in Cai's paper
    o <- order(lfdrs, p)
    lfdrs_sorted <- lfdrs[o]
    fdr_estimate <- cumsum(lfdrs_sorted) / (1:length(p))
    adj_p <- rev(cummin(rev(fdr_estimate)))
    adj_p <- adj_p[order(o)]
    return(adj_p)
}


## helper function for clfdr
lfdr_fit <- function(unadj_p, group, lfdr_estimation = "fdrtool") {
    pvals_list <- split(unadj_p, group)
    if (lfdr_estimation == "fdrtool"){
        lfdr_fun <- function(pv) fdrtool::fdrtool(pv, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
    } else if (lfdr_estimation == "locfdr"){
        if (!requireNamespace("locfdr", quietly=TRUE)){
            stop("locfdr package required for this function to work.")
        }
        lfdr_fun <- function(pv) locfdr::locfdr(qnorm(pv), nulltype=0, plot=0)$fdr
    } else {
        stop("This lfdr estimation method is not available.")
    }
    lfdr_list <- lapply(pvals_list, lfdr_fun)
    lfdrs <- unsplit(lfdr_list, group)

    fit_obj <- data.frame(pvalue=unadj_p, lfdr=lfdrs, group=group)
    fit_obj
}
