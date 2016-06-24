\name{exampleData}
\alias{exampleData}
\docType{data}
\title{GNEPackageMaker example data}
\description{example clinical tab-delimited dataset for GNEPackageMaker package
}
\usage{data(exampleData)}
\format{A tab-delimited table with 4 columns:
 Ids Gender Age Disease
}

\examples{
data(exampleData)

#subset for cancer patients only
cancerPts<-exampleData[which(exampleData$Disease == "Cancer"),]

}
\keyword{datasets}
