\docType{class}
\name{compareClusterResult-class}
\alias{compareClusterResult-class}
\alias{plot,compareClusterResult-method}
\alias{show,compareClusterResult-method}
\alias{summary,compareClusterResult-method}
\title{Class "compareClusterResult"
This class represents the comparison result of gene clusters by GO
categories at specific level or GO enrichment analysis.}
\description{
Class "compareClusterResult" This class represents the
comparison result of gene clusters by GO categories at
specific level or GO enrichment analysis.
}
\section{Slots}{

\describe{
\item{\code{compareClusterResult}}{cluster comparing result}

\item{\code{geneClusters}}{a list of genes}

\item{\code{fun}}{one of groupGO, enrichGO and enrichKEGG}
}}
\author{
Guangchuang Yu \url{http://ygc.name}
}
\seealso{
\code{\linkS4class{groupGOResult}}
\code{\linkS4class{enrichResult}}
\code{\link{compareCluster}}
}
\keyword{classes}

