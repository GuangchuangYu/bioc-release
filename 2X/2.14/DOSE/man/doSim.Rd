\name{doSim}
\alias{doSim}
\title{doSim}
\usage{
  doSim(DOID1, DOID2, measure = "Wang")
}
\arguments{
  \item{DOID1}{DO term vector}

  \item{DOID2}{DO term vector}

  \item{measure}{one of "Wang", "Resnik", "Rel", "Jiang",
  and "Lin".}
}
\value{
  score matrix
}
\description{
  measuring similarities between two DO term vectors.
}
\details{
  provide two DO term vectors, this function will calculate
  their similarities.
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

