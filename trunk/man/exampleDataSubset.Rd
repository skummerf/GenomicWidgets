\name{exampleDataSubset}
\alias{exampleDataSubset}
\docType{data}
\title{subset of data from exampleData.txt}
\description{example subset data object from reading in exampleData.txt}
\usage{data(exampleDataSubset)}
\format{
  A data frame with 2 observations on the following 4 variables.
  \describe{
    \item{\code{Ids}}{a numeric vector}
    \item{\code{Gender}}{a factor with levels \code{female} \code{male}}
    \item{\code{Age}}{a numeric vector}
    \item{\code{Disease}}{a factor with levels \code{Asthma} \code{Cancer}}
  }
}

\examples{
data(exampleDataSubset)
femalePts<-exampleDataSubset[which(exampleDataSubset$Gender == "female"),]
}

\keyword{datasets}
