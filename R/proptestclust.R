##########################################################
## PROPORTION TEST FOR ICS
##
## From 'Proportion test CC2.R'
##########################################################

#' Test of Marginal Proportion for Clustered Data
#'
#' \code{prop.test.clust} can be used for testing the null that the marginal proportion
#' (probability of success) is equal to certain given values in clustered data with potentially informative
#' cluster size.
#' @param x  a vector of binary indicators denoting success/failure of each observation, or a two-dimensional table
#' (or matrix) with 2 columns giving the aggregate counts of failures and successes (respectively) across clusters.
#' @param id a vector which identifies the clusters; ignored if \code{x} is a matrix or table.
#' The length of \code{id} must be the same as the length of \code{x}.
#' @param p the null hypothesized value of the marginal proportion. Must be a single number greater than 0 and less
#' than 1.
#' @param alternative a character string specifying the alternative hypothesis. Must be one of "\code{two.sided}",
#' "\code{greater}", or "\code{less}". You can specify just the initial letter.
#' @param variance character string specifying the method of variance estimation. Must be one of "\code{sand.null}",
#' "\code{sand.est}", "\code{emp}", or "\code{MoM}".
#' @param conf.level confidence level of the returned confidence interval. Must be a single number between
#' 0 and 1.
#' @return A list with class "\code{htest}" containing the following compoments:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{estimate}{the estimated marginal proportion.}
#' \item{null.value}{the value of \code{p} under the null hypothesis.}
#' \item{conf.int}{a confidence interval for the true marginal proportion.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating the test performed and method of construction.}
#' \item{data.name}{a character string giving the name of the data and the total number of clusters.}
#' \item{m}{the number of clusters.}
#' @references
#' Gregg, M., Datta, S., Lorenz, D. (2020) Variance Estimation in Tests of Clustered Categorical Data
#' with Informative Cluster Size.
#' \emph{Submitted}, \bold{xx}, xx--xx.
#'
#' @details If \code{p} is not given, the null tested is that the underlying marginal probablity of
#' success is .5.
#'
#' The \code{variance} argument allows the user to specify the method of variance estimation, selecting from
#' the sandwich estimate evaluated at the null hypothesis (\code{sand.null}), the sandwich estimate evaluated at the
#' cluster-weighted proportion (\code{sand.est}), the emperical estimate (\code{emp}), or the method of moments
#' estimate (\code{MoM}).
#'
#' @examples
#' data(icsdat)
#' prop.test.clust(icsdat$cat1, icsdat$id)
#' prop.test.clust(table(icsdat$id, icsdat$cat1), variance="emp")
#'
#' @export


prop.test.clust <- function(x, id, p = NULL, alternative = c("two.sided", "less", "greater"),
                            variance=c("sand.null", "sand.est", "emp", "MoM"), conf.level = 0.95) {
  alternative <- match.arg(alternative)
  variance <- match.arg(variance)
  DNAME <- deparse(substitute(x))

  if (is.table(x)) {
    if (ncol(x) != 2L)
      stop("'x' must have 2 columns")
    x <- x[complete.cases(x),]
    if (!(all(rowSums(x)>0)))
      stop("all clusters must have counts > 0")
    if (!(all(x)>=0))
      stop("elements of x must be nonnegative ")
    m <- nrow(x)
    xbar <- x[,2]/rowSums(x)
    xtot <- sum(x[,2])
    ntot <- sum(x)
  }
  else {
    if ((l <- length(x)) != length(id))
      stop("'x' and 'id' must have the same length")
    OK <- complete.cases(x, id)
    x <- x[OK]
    x <- 1*x
    id <- id[OK]
    id <- factor(as.factor(id), levels=unique(as.factor(id)))
    if ((k <- length(x)) < 1L)
      stop("not enough data")
    if (!all(x%in%c(0,1)))
      stop("elements of 'x' must be binary")
    m <- length(unique(id))
    xbar <- aggregate(x, list(id), mean)[,2]
    xtot <- sum(x)
    ntot <- length(x)
  }

  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if ((!is.null(p)) && (p <= 0) || (!is.null(p)) && (p >= 1))
    stop("elements of 'p' must be in (0,1)")
  if (is.null(p))
    p <- 0.5


  DNAME <- paste0(paste0(DNAME, ", M = "), as.character(m))
  alpha <- 1-conf.level
  p.hat <- mean(xbar)

  ##
  ## Variance
  H <- switch(variance, sand.null = (-p.hat/p^2 - (1-p.hat)/(1-p)^2), sand.est = (-1/p.hat)-(1/(1-p.hat)),
              MoM = 1, emp = 1)

  V <- switch(variance, sand.null = mean(((xbar-p)/(p*(1-p)))^2),
              sand.est = mean(((xbar-p.hat)/(p.hat*(1-p.hat)))^2), MoM = 1, emp = 1)

  var.hat <- switch(variance, sand.null = (1/m)*(H^-1)*V*(H^-1), sand.est = (1/m)*(H^-1)*V*(H^-1),
                    MoM = (1/m)*mean((xbar-p)^2), emp = (1/m)*var(xbar))


  z <- (p.hat - p)/(sqrt(var.hat))

  pval <- switch(alternative, less = pnorm(z), greater = pnorm(z, lower.tail = FALSE),
                 two.sided = 2*(1-pnorm(abs(z))))

  cint <- switch(alternative, less = c(0, min(p.hat+qnorm(conf.level)*sqrt(var.hat),1)),
                 greater = c(max(p.hat - qnorm(conf.level)*sqrt(var.hat),0), 1),
                 two.sided = c(max(p.hat - qnorm(1-alpha/2)*sqrt(var.hat),0),
                               min(p.hat + qnorm(1-alpha/2)*sqrt(var.hat),1)))

  METHOD <- paste0("Cluster-weighted proportion test with variance est: ", variance)

  names(m) <- "M"
  names(p.hat) <- "Cluster-weighted proportion"
  names(p) <- "p"
  names(z) <- "z"
  attr(cint, "conf.level") <- conf.level

  RVAL <- list(statistic = z,
               p.value = pval, estimate = p.hat, null.value = p,
               conf.int = cint, alternative = alternative, method = METHOD,
               data.name = DNAME, m=m)
  class(RVAL) <- "htest"
  if (m < 30)
    warning('Number of clusters < 30. Normal approximation may be incorrect')
  return(RVAL)
}

