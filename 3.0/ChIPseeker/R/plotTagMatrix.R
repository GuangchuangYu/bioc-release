##' plot the profile of peaks
##'
##' 
##' @title plotAvgProf 
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param xlim xlim
##' @param xlab x label
##' @param ylab y label
##' @return ggplot object
##' @export
##' @author G Yu
plotAvgProf <- function(tagMatrix, xlim,
                        xlab="Genomic Region (5'->3')",
                        ylab = "Read Count Frequency") {
    p <- plotAvgProf.internal(tagMatrix, xlim=xlim,
                              xlab=xlab, ylab=ylab)
    
    return(p)
}

##' plot the profile of peaks that align to flank sequences of TSS 
##'
##' 
##' @title plotAvgProf
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param xlab xlab
##' @param ylab ylab
##' @param verbose print message or not
##' @return ggplot object
##' @export
##' @author G Yu
plotAvgProf2 <- function(peak, weightCol=NULL, TxDb=NULL,
                        upstream=1000, downstream=1000,
                        xlab="Genomic Region (5'->3')",
                        ylab="Read Count Frequency",
                        verbose=TRUE) {
    
    if (verbose) {
        cat(">> preparing promoter regions...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    promoter <- getPromoters(TxDb=TxDb,
                             upstream=upstream,
                             downstream=downstream)
    
    if (verbose) {
        cat(">> preparing tag matrix...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    if ( is(peak, "list") ) {
        tagMatrix <- lapply(peak, getTagMatrix,
                            weightCol=weightCol, windows=promoter)
    } else {
        tagMatrix <- getTagMatrix(peak, weightCol, promoter)
    }
    
    if (verbose) {
        cat(">> plotting figure...\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    p <- plotAvgProf.internal(tagMatrix, xlim=c(-upstream, downstream),
                               xlab=xlab, ylab=ylab)
    return(p)
}

##' plot the heatmap of tagMatrix
##'
##' 
##' @title tagHeatmap
##' @param tagMatrix tagMatrix or a list of tagMatrix
##' @param xlim xlim
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param color color
##' @return figure
##' @export
##' @author G Yu
tagHeatmap <- function(tagMatrix, xlim, xlab="", ylab="", title=NULL, color="red") {
    listFlag <- FALSE
    if (is(tagMatrix, "list")) {
        listFlag <- TRUE
    }
    peakHeatmap.internal2(tagMatrix, xlim, listFlag, color, xlab, ylab, title)
}

##' plot the heatmap of peaks align to flank sequences of TSS
##'
##' 
##' @title peakHeatmap
##' @param peak peak file or GRanges object
##' @param weightCol column name of weight
##' @param TxDb TxDb object
##' @param upstream upstream position
##' @param downstream downstream position
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @param color color
##' @param verbose print message or not
##' @return figure
##' @export
##' @author G Yu
peakHeatmap <- function(peak, weightCol=NULL, TxDb=NULL,
                            upstream=1000, downstream=1000,
                            xlab="", ylab="", title=NULL,
                            color=NULL, verbose=TRUE) {
    listFlag <- FALSE
    if ( is(peak, "list") ) {
        listFlag <- TRUE
        if (is.null(names(peak)))
            stop("peak should be a peak file or a name list of peak files...")
    }
        
    if (verbose) {
        cat(">> preparing promoter regions...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    promoter <- getPromoters(TxDb=TxDb,
                             upstream=upstream, downstream=downstream)

    if (verbose) {
        cat(">> preparing tag matrix...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    if (listFlag) {
        tagMatrix <- lapply(peak, getTagMatrix, weightCol=weightCol, windows=promoter)
    } else {
        tagMatrix <- getTagMatrix(peak, weightCol, promoter)
    }

    if (verbose) {
        cat(">> generating figure...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
   
    xlim=c(-upstream, downstream)
    
    peakHeatmap.internal2(tagMatrix, xlim, listFlag, color, xlab, ylab, title)

    if (verbose) {
        cat(">> done...\t\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }
    invisible(tagMatrix)
}

peakHeatmap.internal2 <- function(tagMatrix, xlim, listFlag, color, xlab, ylab, title) {
    if ( is.null(xlab) || is.na(xlab))
        xlab <- ""
    if ( is.null(ylab) || is.na(ylab))
        ylab <- ""

    if (listFlag) {
        nc <- length(tagMatrix)
        if ( is.null(color) || is.na(color) ) {
            cols <- getCols(nc)
        } else if (length(color) != nc) {
            cols <- rep(color[1], nc)
        } else {
            cols <- color
        }
        
        if (is.null(title) || is.na(title))
            title <- names(tagMatrix)
        if (length(xlab) != nc) {
            xlab <- rep(xlab[1], nc)
        }
        if (length(ylab) != nc) {
            ylab <- rep(ylab[1], nc)
        }
        if (length(title) != nc) {
            title <- rep(title[1], nc)
        }
        par(mfrow=c(1, nc))
        for (i in 1:nc) {
            peakHeatmap.internal(tagMatrix[[i]], xlim, cols[i], xlab[i], ylab[i], title[i])
        }
    } else {
        if (is.null(color) || is.na(color))
            color <- "red"
        if (is.null(title) || is.na(title))
            title <- ""
        peakHeatmap.internal(tagMatrix, xlim, color, xlab, ylab, title)
    }
}


##' @importFrom grDevices colorRampPalette
##' @importFrom BiocGenerics image
peakHeatmap.internal <- function(tagMatrix, xlim=NULL, color="red", xlab="", ylab="", title="") {
    tagMatrix <- t(apply(tagMatrix, 1, function(x) x/max(x)))
    ii <- order(rowSums(tagMatrix))
    tagMatrix <- tagMatrix[ii,]
    cols <- colorRampPalette(c("white",color))(200)    
    if (is.null(xlim)) {
        xlim <- 1:ncol(tagMatrix)
    } else if (length(xlim) == 2) {
        xlim <- seq(xlim[1], xlim[2])
    } 
    image(x=xlim, y=1:nrow(tagMatrix),z=t(tagMatrix),useRaster=TRUE, col=cols, yaxt="n", ylab="", xlab=xlab, main=title)
}


##' @importFrom plyr ldply
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
plotAvgProf.internal <- function(tagMatrix,
                                  xlim=c(-3000,3000),
                                  xlab="Genomic Region (5'->3')",
                                  ylab="Read Count Frequency") {

    listFlag <- FALSE
    if (is(tagMatrix, "list")) {
        if ( is.null(names(tagMatrix)) ) {
            stop("tagMatrix should be a named list...")
        }
        listFlag <- TRUE
    } 
    
    if ( listFlag ) {
        if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix[[1]]) ) {
            stop("please specify appropreate xcoordinations...")
        }
    } else {
        if ( (xlim[2]-xlim[1]+1) != ncol(tagMatrix) ) {
            stop("please specify appropreate xcoordinations...")
        }
    }

    pos <- value <- .id <- NULL
    
    if ( listFlag ) {
        tagCount <- lapply(tagMatrix, getTagCount, xlim=xlim)
        tagCount <- ldply(tagCount)
        p <- ggplot(tagCount, aes(pos, value, group=.id, color=.id))
    } else {
        tagCount <- getTagCount(tagMatrix, xlim=xlim)
        p <- ggplot(tagCount, aes(pos, value))
    }

    p <- p + geom_line()
    if ( 0 > xlim[1] && 0 < xlim[2] ) {
        p <- p + geom_vline(xintercept=0,
                            linetype="longdash")
        p <- p + scale_x_continuous(breaks=c(xlim[1], floor(xlim[1]/2),
                                       0,
                                       floor(xlim[2]/2), xlim[2]),
                                   labels=c(xlim[1], floor(xlim[1]/2),
                                       "TSS",
                                       floor(xlim[2]/2), xlim[2]))
    }

    if (listFlag) {
        cols <- getCols(length(tagMatrix))
        p <- p + scale_color_manual(values=cols)
    }
    p <- p+xlab(xlab)+ylab(ylab)
    p <- p + theme_bw() + theme(legend.title=element_blank())
    return(p)
}
