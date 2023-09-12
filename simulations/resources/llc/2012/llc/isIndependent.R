isIndependent <-function( C, x, y, given=NULL, N ) {
  # Returns a p-value of an independence test.
  #INPUT:
  # C         - Covariance matrix of the data.
  # x,y       - Indexes of variables in question.
  # given     - A vector of variable indexes that are conditioned on.
  # N         - The number of samples.
  #OUTPUT:
  # p-value of the test

  #calculating first the partial correlation
  p <- pcor( c(x,y,given), C)

  #returning a p-value of this test 
  if ( is.infinite(N) ) {
    N = 1e15
  }
  pcor.test(p, length(given), N )$pvalue
}

# From https://github.com/cran/ggm/blob/master/R/functions.R (Lines 2784-2848)
#' Partial correlation
#' 
#' Computes the partial correlation between two variables given a set of other
#' variables.
#' 
#' 
#' @param u a vector of integers of length > 1. The first two integers are the
#' indices of variables the correlation of which must be computed. The rest of
#' the vector is the conditioning set.
#' @param S a symmetric positive definite matrix, a sample covariance matrix.
#' @return a scalar, the partial correlation matrix between variables
#' \code{u[1]} and \code{u[2]} given \code{u[-c(1,2)]}.
#' @author Giovanni M. Marchetti
#' @seealso \code{\link{cor}}, \code{\link{parcor}}, \code{\link{correlations}}
#' @keywords models multivariate
#' @examples
#' 
#' data(marks)
#' ## The correlation between vectors and algebra given analysis and statistics
#'  pcor(c("vectors", "algebra", "analysis", "statistics"), var(marks))
#' ## The same
#' pcor(c(2,3,4,5), var(marks))
#' ## The correlation between vectors and algebra given statistics
#'  pcor(c("vectors", "algebra", "statistics"), var(marks))
#' ## The marginal correlation between analysis and statistics 
#' pcor(c("analysis","statistics"), var(marks))
#' 
"pcor" <-
function (u, S) 
{
### Partial correlation between u[1:2], given th rest of u. S: cov matrix.
  k <- solve(S[u,u])
  -k[1,2]/sqrt(k[1,1]*k[2,2])
}



#' Test for zero partial association
#' 
#' Test for conditional independence between two variables, given the other
#' ones, assuming a multivariate normal distribution.
#' 
#' 
#' @param r a partial correlation coefficient, computed by \code{\link{pcor}}.
#' @param q the number of variables in the conditioning set.
#' @param n integer > 0, the sample size.
#' @return \item{tval}{The Student's t-test statistic.} \item{df}{The degrees
#' of freedom} \item{pvalue}{The P-value, assuming a two-sided alternative.}
#' @author Giovanni M. Marchetti
#' @seealso \code{\link{pcor}}, \code{\link{shipley.test}}
#' @keywords htest multivariate
#' @examples
#' 
#' ## Are 2,3 independent given 1?
#' data(marks)
#' pcor.test(pcor(c(2,3,1), var(marks)), 1, n=88)
#' 
"pcor.test" <-
function(r, q, n){
                df = n - 2 - q
                tval <- r * sqrt(df)/sqrt(1-r*r)
                pv <- 2 * pt(-abs(tval), df)
  list(tval = tval, df = df, pvalue = pv)

}
