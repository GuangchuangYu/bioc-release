##' convert a list of gene IDs to igraph0 object.
##'
##'
##' @title convert gene IDs to igraph0 object
##' @param inputList a list of gene IDs
##' @return a igraph0 object.
##' @importFrom igraph0 graph.data.frame
##' @author Guangchuang Yu \url{http://ygc.name}
list2graph <- function(inputList) {
    x <- data.frame()
    for (i in 1:length(inputList)) {
        x=rbind(x,
        data.frame(categoryID=rep(names(inputList[i]),
                   length(inputList[[i]])),
                   Gene=inputList[[i]]))
    }
    g <- graph.data.frame(x, directed=F)
    return(g)
}


##' plot function of gene Concept Net.
##'
##'
##' @title plot gene net by categories
##' @param inputList a list of gene IDs
##' @param categorySize setting category size
##' @param showCategory number of categories to plot
##' @param pvalue pvalue
##' @param logFC  log fold Change
##' @param output output type
##' @return plotted igraph0 object.
##' @importFrom scales cscale
##' @importFrom scales seq_gradient_pal
##' @importFrom igraph0 tkplot
##' @importFrom igraph0 plot.igraph
##' @importFrom igraph0 V
##' @importFrom igraph0 "V<-"
##' @importFrom igraph0 E
##' @importFrom igraph0 "E<-"
##' @importFrom igraph0 degree
##' @importFrom igraph0 layout.fruchterman.reingold
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
cnetplot <- function(inputList, categorySize="geneNum",
                     showCategory=5, pvalue=NULL,
                     logFC=NULL, output="fixed") {

    inputList <- inputList[1:showCategory]
    pvalue <- pvalue[1:showCategory]

    ## generate graph object
    g=list2graph(inputList)

    ## setting some attributes
    V(g)$size <- 8
    V(g)$color <- "#B3B3B3"
    V(g)$label <- V(g)$name

    E(g)$width=2
    E(g)$color <- "#8DA0CB"

    ## attributes of Category Node
    lengthOfCategory <- length(inputList)
    V(g)[0:(lengthOfCategory-1)]$size=30  ## setting by default.
    V(g)[0:(lengthOfCategory-1)]$color= "#E5C494"

    if(is.numeric(categorySize)) {
        V(g)[0:(lengthOfCategory-1)]$size=categorySize
    } else {
        if (categorySize == "geneNum") {
            n <- degree(g)[1:lengthOfCategory]
            V(g)[0:(lengthOfCategory-1)]$size <- n/sum(n) * 100
        }
        if (categorySize == "pvalue") {
            if (is.null(pvalue)) {
                stop("pvalue must not be null...")
            }
            pScore<- -log10(pvalue)
            V(g)[0:(lengthOfCategory-1)]$size <- pScore/sum(pScore) * 100
        }
    }


    if ((!is.null(logFC)) &
        (all(unique(unlist(inputList)) %in% names(logFC)))) {

        logFC.down <- logFC[logFC < 0]
        logFC.up <- logFC[logFC >=0]

        col.down <- cscale(logFC.down, seq_gradient_pal("darkgreen", "green"))
        col.up <- cscale(logFC.up, seq_gradient_pal("red", "darkred"))
        col <- c(col.down, col.up)


        gn <- V(g)[lengthOfCategory:length(V(g))]$label
        V(g)[lengthOfCategory:length(V(g))]$color = col[gn]
    }


    if (output == "fixed"){
        plot.igraph(g,
                    vertex.label.font=2,
                    vertex.label.color='#666666',
                    vertex.label.cex=1.5,
                    vertex.frame.color=V(g)$color,
                    layout=layout.fruchterman.reingold)
    } else {
        tkplot(g,
               vertex.label.font = 2,
               vertex.label.color = '#666666',
               vertex.label.cex=1.5,
               vertex.frame.color=V(g)$color,
               layout=layout.fruchterman.reingold)
    }

}

cnetplot.enrichResult <- function(x,
                                  showCategory=5,
                                  categorySize="geneNum",
                                  pvalue=NULL,
                                  logFC=NULL,
                                  output="fixed") {
    res <- summary(x)
    gc <- x@geneInCategory

    if ("pvalue" %in% names(res)) {
        y <- res[res$ID %in% names(gc),
                 c("ID", "Description", "pvalue")]
    } else {
        y <- res[res$ID %in% names(gc),
                 c("ID", "Description")]

    }

    gc <- gc[as.character(y$ID)]
    names(gc) <- y$Description

    if (showCategory > length(gc)) {
        showCategory = length(gc)
    }

    if (categorySize == "pvalue") {
        pvalue <- y$pvalue
    } else {
        pvalue <- NULL
    }

    readable <- x@readable
    organism <- x@organism
    if (readable & (!is.null(logFC) ) ){
        gn <- EXTID2NAME(names(logFC),organism)
        names(logFC) <- gn
    }

    cnetplot(inputList=gc,
             showCategory=showCategory,
             categorySize=categorySize,
             pvalue=pvalue,
             logFC=logFC,
             output=output)
}
