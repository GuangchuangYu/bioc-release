\name{enrichPathway}
\alias{enrichPathway}
\title{Pathway Enrichment Analysis of a gene set.
Given a vector of genes, this function will return the enriched pathways
with FDR control.}
\usage{
  enrichPathway(gene, organism = "human",
    pvalueCutoff = 0.05, pAdjustMethod = "BH",
    qvalueCutoff = 0.2, universe, minGSSize = 5,
    readable = FALSE)
}
\arguments{
  \item{gene}{a vector of entrez gene id.}

  \item{organism}{one of "human", "rat", "mouse",
  "celegans", "zebrafish", "fly".}

  \item{pvalueCutoff}{Cutoff value of pvalue.}

  \item{pAdjustMethod}{one of "holm", "hochberg", "hommel",
  "bonferroni", "BH", "BY", "fdr", "none"}

  \item{qvalueCutoff}{Cutoff value of qvalue}

  \item{universe}{background genes}

  \item{minGSSize}{minimal size of genes annotated by
  Ontology term for testing.}

  \item{readable}{whether mapping gene ID to gene Name}
}
\value{
  A \code{enrichResult} instance.
}
\description{
  Pathway Enrichment Analysis of a gene set. Given a vector
  of genes, this function will return the enriched pathways
  with FDR control.
}
\examples{
gene <- c("11171", "8243", "112464", "2194",
				"9318", "79026", "1654", "65003",
				"6240", "3476", "6238", "3836",
				"4176", "1017", "249")
	yy = enrichPathway(gene, pvalueCutoff=0.05)
	head(summary(yy))
	#plot(yy)
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\seealso{
  \code{\link{enrichResult-class}}
}
\keyword{manip}

