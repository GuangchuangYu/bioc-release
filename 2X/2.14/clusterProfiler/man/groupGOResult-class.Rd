\docType{class}
\name{groupGOResult-class}
\alias{groupGOResult-class}
\alias{show,groupGOResult-method}
\title{Class "groupGOResult"
This class represents the result of functional Profiles of a set of gene at
specific GO level.}
\description{
Class "groupGOResult" This class represents the result of
functional Profiles of a set of gene at specific GO level.
}
\section{Slots}{

\describe{
\item{\code{result}}{GO classification result}

\item{\code{ontology}}{Ontology}

\item{\code{level}}{GO level}

\item{\code{organism}}{one of "human", "mouse" and "yeast"}

\item{\code{gene}}{Gene IDs}

\item{\code{geneInCategory}}{gene and category association}

\item{\code{readable}}{logical flag of gene ID in symbol or not.}
}}
\author{
Guangchuang Yu \url{http://ygc.name}
}
\seealso{
\code{\linkS4class{compareClusterResult}}
\code{\link{compareCluster}} \code{\link{groupGO}}
}
\keyword{classes}

