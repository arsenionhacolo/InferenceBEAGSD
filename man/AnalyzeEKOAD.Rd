% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalyzeEKOAD.r
\name{AnalyzeEKOAD}
\alias{AnalyzeEKOAD}
\title{Analyze simulated adaptive trials.}
\usage{
AnalyzeEKOAD(replicates = NULL, basedir = NULL)
}
\arguments{
\item{replicates}{Number of simulated trials to be analysed. If \code{NULL} (default),
all trials found in \code{./basedir/SimulatedTrials} are analysed.}

\item{basedir}{The base directory containing the sub-directory \code{SimulatedTrials}
with the simulated trials. If \code{NULL} (default), the current working directory
is uded.}
}
\value{
A dataframe with the results. A copy is saved in the file \code{Results.csv}
in the \code{basedir}.
}
\description{
\code{AnalyzeEKOAD} performs inference on trials simulated by the function
\code{\link{SimulateEKOAD}} using the methods proposed by Nhacolo and Brannath (2018)
and naive maximum likelihood.
}
\details{
Overall p-values, point estimates and confidence intervals are calculated.
}
\references{
Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
\emph{Stat Methods Med Res}, 2018.
}
\seealso{
\code{\link{SimulateEKOAD}, \code{\link{mue1}}, \code{\link{mue2}}, \code{\link{mue2v2}}, \code{\link{mue3}}}.
}
\author{
Arsenio Nhacolo
}
