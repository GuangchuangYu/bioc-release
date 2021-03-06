\name{scaleNodeColor}
\alias{scaleNodeColor}
\title{scaleNodeColor}
\usage{
  scaleNodeColor(g, foldChange, node.idx = NULL,
    DE.foldChange)
}
\arguments{
  \item{g}{igraph object}

  \item{foldChange}{fold Change}

  \item{node.idx}{index of node to color}

  \item{DE.foldChange}{logical}
}
\value{
  igraph object
}
\description{
  scale color nodes
}
\details{
  color nodes based on fold change of expression
}
\author{
  Yu Guangchuang
}

