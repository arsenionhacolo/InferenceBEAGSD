% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvaluek.r
\name{pvaluek}
\alias{pvaluek}
\title{P-value}
\usage{
pvaluek(s, n1, r1, n, p0)
}
\arguments{
\item{s}{Total number of successes.}

\item{n1}{Stage 1 sample size.}

\item{r1}{Stage 1 critical value (trial is stopped at stage 1 if the number of successes
is at most \code{r1}).}

\item{n}{Total sample size.}

\item{p0}{Response rate under the null hypothesis.}
}
\value{
P-value.
}
\description{
Calculates p-value for Simon-like designs.
}
\details{
It is based on the definition of p-value by \emph{Koyama and Chen (2008)}.
}
\examples{
pvaluek(21, 19, 4, 54, 0.2)
}
\references{
Koyama, T. and Chen, H. Proper inference from Simon's two-stage designs.
\emph{Stat Med}, 2008, 27, 3145-3154.
}
\seealso{
\code{\link{pquantile}} and \code{\link{pk}}.
}
\author{
Arsenio Nhacolo
}
