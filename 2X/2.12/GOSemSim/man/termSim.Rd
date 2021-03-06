\name{termSim}
\alias{termSim}
\title{termSim}
\usage{
  termSim(t1, t2, method = "Wang", organism = "human", ont)
}
\arguments{
  \item{t1}{term vector}

  \item{t2}{term vector}

  \item{method}{one of "Wang", "Resnik", "Rel", "Jiang",
  and "Lin".}

  \item{organism}{only "human" supported}

  \item{ont}{ontology}
}
\value{
  score matrix
}
\description{
  measuring similarities between two term vectors.
}
\details{
  provide two term vectors, this function will calculate
  their similarities.
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

