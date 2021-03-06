% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/plotAnno.R
\name{plotAnnoBar.data.frame}
\alias{plotAnnoBar.data.frame}
\title{plotAnnoBar.data.frame}
\usage{
plotAnnoBar.data.frame(anno.df, xlab = "", ylab = "Percentage(\%)",
  title = "Feature Distribution", categoryColumn)
}
\arguments{
\item{anno.df}{annotation stats}

\item{xlab}{xlab}

\item{ylab}{ylab}

\item{title}{plot title}

\item{categoryColumn}{category column}
}
\value{
bar plot that summarize genomic features of peaks
}
\description{
plot feature distribution based on their chromosome region
}
\details{
plot chromosome region features
}
\author{
Guangchuang Yu \url{http://ygc.name}
}
\seealso{
\code{\link{annotatePeak}} \code{\link{plotAnnoPie}}
}

