##' Compare gene clusters functional profile
##' Given a list of gene set, this function will compute profiles of each gene
##' cluster.
##'
##'
##' @param geneClusters a list of entrez gene id.
##' @param fun One of groupGO and enrichGO.
##' @param ...  Other arguments.
##' @return A \code{clusterProfResult} instance.
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{compareClusterResult-class}}, \code{\link{groupGO}}
##'   \code{\link{enrichGO}}
##' @keywords manip
##' @examples
##'
##' 	data(gcSample)
##' 	xx <- compareCluster(gcSample, fun=enrichKEGG, organism="human", pvalueCutoff=0.05)
##' 	#summary(xx)
##' 	#plot(xx, type="dot", caption="KEGG Enrichment Comparison")
##'
compareCluster <- function(geneClusters, fun=enrichGO, ...) {
    clProf <- llply(geneClusters,
                    .fun=function(i) {
			x=fun(i, ...)
			summary(x)
                    }
                    )

    clProf.df <- ldply(clProf, rbind)
    ##colnames(clProf.df)[1] <- "Cluster"
    clProf.df <- rename(clProf.df, c(.id="Cluster"))
    new("compareClusterResult",
        compareClusterResult = clProf.df,
        geneClusters = geneClusters,
        fun = fun
	)
}


##' Class "compareClusterResult"
##' This class represents the comparison result of gene clusters by GO
##' categories at specific level or GO enrichment analysis.
##'
##'
##' @name compareClusterResult-class
##' @aliases compareClusterResult-class show,compareClusterResult-method
##'   summary,compareClusterResult-method plot,compareClusterResult-method
##' @docType class
##' @slot compareClusterResult cluster comparing result
##' @slot geneClusters a list of genes
##' @slot fun one of groupGO, enrichGO and enrichKEGG
##' @exportClass compareClusterResult
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\linkS4class{groupGOResult}}
##'   \code{\linkS4class{enrichGOResult}} \code{\link{compareCluster}}
##' @keywords classes
setClass("compareClusterResult",
         representation = representation(
         compareClusterResult = "data.frame",
         geneClusters = "list",
         fun = "function"
         )
         )

##' show method for \code{compareClusterResult} instance
##'
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @param object A \code{compareClusterResult} instance.
##' @return message
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="compareClusterResult"),
          function (object){
              geneClusterLen <- length(object@geneClusters)
              #fun <- object@fun
              #fun <- as.character(substitute(fun))
              #if (fun == "enrichKEGG") {
  #                analysis <- "KEGG Enrichment Analysis"
  #            } else if (fun == "groupGO") {
  #                analysis <- "GO Profiling Analysis"
  #           } else if (fun == "enrichGO") {
  #                analysis <- "GO Enrichment Analysis"
   #           } else if (fun == "enrichDO") {
    #              analysis <- "DO Enrichment Analysis"		 
#	      } else {
#		analysis <- "User specify Analysis"
#	      }
#              cat ("Compare", geneClusterLen, "gene clusters using", analysis, "\n")
		cat ("Result of Comparing", geneClusterLen, "gene clusters", "\n")
          }
          )

##' summary method for \code{compareClusterResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{compareClusterResult} instance.
##' @return A data frame
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("summary", signature(object="compareClusterResult"),
          function(object) {
              return(object@compareClusterResult)
          }
          )

##' plot method for \code{compareClusterResult} instance
##'
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##'
##' @title plot method
##' @param x A \code{compareClusterResult} instance.
##' @param type one of "dot" and "bar".
##' @param title graph title
##' @param font.size graph font size
##' @param limit numeric parameter, restrict the top categories for plotting.
##' @param by one of "percentage" and "count"
##' @return ggplot object
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(x="compareClusterResult"),
          function(x, type="dot", title="", font.size=12, limit=5, by="percentage") {
              clProf.df <- summary(x)

              ## get top 5 (default) categories of each gene cluster.
              if (is.null(limit)) {
                  result <- clProf.df
              } else {
                  Cluster <- NULL # to satisfy codetools
                  result <- ddply(.data = clProf.df,
                                  .variables = .(Cluster),
                                  .fun = function(df, N) {
                                      if (length(df$Count) > N) {
                                          idx <- order(df$Count, decreasing=T)[1:N]
                                          return(df[idx,])
                                      } else {
                                          return(df)
                                      }
                                  },
                                  N=limit
                                  )
              }

              ## remove zero count
              result$Description <- as.character(result$Description) ## un-factor
              GOlevel <- result[,c(2,3)] ## GO ID and Term
              GOlevel <- unique(GOlevel)

              result <- result[result$Count != 0, ]
              result$Description <- factor(result$Description, levels=rev(GOlevel[,2]))


              if (by=="percentage") {
                  Description <- Count <- NULL # to satisfy codetools
                  result <- ddply(result, .(Description), transform, Percentage = Count/sum(Count), Total = sum(Count))

                  ## label GO Description with gene counts.
                  x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
                  y <- sapply(x[,3], paste, ")", sep="")
                  result$Description <- y

                  ## restore the original order of GO Description
                  xx <- result[,c(2,3)]
                  xx <- unique(xx)
                  rownames(xx) <- xx[,1]
                  Termlevel <- xx[as.character(GOlevel[,1]),2]

                  ##drop the *Total* column
                  result <- result[, colnames(result) != "Total"]

                  result$Description <- factor(result$Description, levels=rev(Termlevel))

              } else if (by == "count") {

              } else {

              }
              p <- plotting.clusterProfile(result, type, by, title, font.size)
			  return(p)
          }
          )
