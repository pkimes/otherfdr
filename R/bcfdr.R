#' Barber-Candès False Discovery Rate
#'
#' @description
#' Simple implementation of the false discovery rate controlling procedure of
#' Barber and Candès (2015). The implementation is based on code available on
#' \href{https://github.com/lihualei71/AdaPT/blob/38a7d0a/R/other_methods.R}{GitHub}
#' by the authors of the \pkg{adaptMT} package as part of their analysis repository.
#' 
#' The function has been modified to return q-values for each p-value rather than
#' simply return the indices of rejected hypotheses at a pre-specified nominal FDR
#' threshold. The q-values are computed as the smallest nominal FDR threshold at
#' which each p-value would be rejected.
#' 
#' @param p numeric vector of unadjusted p-values.
#'
#' @return
#' Numeric vector of q-values of equal length and order
#' as the to input vector of p-values.
#'
#' @references
#' References for the knockoff procedure for false discovery rate control (underlying theory):
#' 
#' Barber, R.F., Candès, E.J. (2015) Controlling the false discovery rate via knockoffs. The Annals of Statistics, 43(5):2055-2085. \url{https://doi.org/10.1214/15-AOS1337}
#'
#' Arias-Castro, E., Chen, S. (2017) Distribution-free multiple testing. Electronic Journal of Statistics, 11(1):1983-2001. \url{https://doi.org/10.1214/17-EJS1277}
#'
#' Reference for the original BC procedure function (underlying code):
#'
#' Lihua, L. (2018) AdaPT paper repoistory. \url{https://github.com/lihualei71/AdaPT/}
#'
#' @examples
#' ## generate example p-values
#' pv <- runif(100)
#'
#' ## q-values using the BC procedure
#' qv <- bcfdr(p = pv)
#' 
#' @seealso \code{\link[knockoff]{knockoff}}
#' @export
#' @author Patrick Kimes
bcfdr <- function(p) {
    sorted.mask.p <- sort(pmin(p, 1 - p))
    ## fdphat only changes at mask point cutoffs 
    fdphat <- sapply(sorted.mask.p, function(thresh) {
        (1 + sum(p >= 1 - thresh)) / max(1, sum(p <= thresh))
    })
    ## determine min mask/cutoff where each hypothesis is rejected
    sapply(p, function(x) { min(c(fdphat[sorted.mask.p >= x], Inf)) })
}
