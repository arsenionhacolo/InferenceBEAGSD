% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pu.r
\name{pu}
\alias{pu}
\title{UMVUE}
\usage{
pu(s, n1, r1, n)
}
\arguments{
\item{s}{Total number of successes.}

\item{n1}{Stage 1 sample size.}

\item{r1}{Stage 1 critical value (trial is stopped at stage 1 if the number of successes
is at most \code{r1}).}

\item{n}{Total sample size.}
}
\value{
Estimate of the response rate.
}
\description{
Calculates the uniformly minimum variance unbiased estimator (UMVUE) of the true
response probability.
}
\details{
The UMVUE is based on approach by \emph{Grishick et al. (1946)}. It was first
considered by \emph{Chang et al. (1989)} and further studied by \emph{Jung et al. (2004)}.
}
\examples{
pu(21, 19, 4, 54)
}
\references{
Jung, S.-H. and Kim, K. M. On the estimation of the binomial probability in
multistage clinical trials. \emph{Stat Med}, 2004, 23, 881-896.
}
\seealso{
\code{\link{pm}}, \code{\link{pg}}, \code{\link{pp}} and
\code{\link{pk}}.
}
\author{
Arsenio Nhacolo
}
