#' Boca-Leek False Discovery Rate
#'
#' @description
#' A simple wrapper function to the false discovery rate controlling procedure of
#' Boca and Leek (2018) which returns adjusted p-values. The heavy lifting of
#' modeling the relationship between covariates and the proportion of null hypotheses
#' is handled by the \code{\link[swfdr]{lm_pi0}} function in the \pkg{swfdr} package.
#' 
#' @param p numeric vector of unadjusted p-values.
#' @param x matrix of covariates with rows corresponding to p-values \code{p};
#'        alternatively, numeric vector for a single covariate of same length
#'        as \code{p}. 
#' @param ... additional parameters to pass to \code{\link[swfdr]{lm_pi0}}.
#' 
#' @return
#' Numeric vector of adjusted p-values of equal length and order
#' as the to input vector of p-values.
#'
#' @references
#' Reference for the false discovery rate controlling procedure (underlying theory):
#' 
#' Boca, S.M., Leek, J.T. (2018) A direct approach to estimating false discovery rates conditional on covariates. bioRxiv. \url{https://doi.org/10.1101/035675}
#' 
#' Reference for the primary null proportion estimating function (underlying code):
#' 
#' Leek, J.T., Jager, L., Boca, S.M. (2018) swfdr: Science-wise false discovery rate and proportion of true null hypotheses estimation. R package version 1.6.0. \url{https://doi.org/doi:10.18129/B9.bioc.swfdr}
#'
#' @examples
#' ## generate example p-values and covariate values
#' pv <- runif(100)
#' x <- sample(100)
#'
#' ## adjust p-values, grouping by covariate bin
#' adjpv <- blfdr(p = pv, x = x)
#'
#' @seealso \code{\link[swfdr]{lm_pi0}}
#' @importFrom swfdr lm_pi0
#' @export
#' @author Patrick Kimes
blfdr <- function(p, x, ...) {
    x <- swfdr::lm_pi0(pValues = p, X = x, ...)
    adj_p <- x$pi0 * p.adjust(pval, method = "BH")
    return(adj_p)
}
