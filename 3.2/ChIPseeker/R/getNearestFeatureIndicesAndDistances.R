##' get index of features that closest to peak and calculate distance
##'
##' 
##' @title getNearestFeatureIndicesAndDistances 
##' @param peaks peak in GRanges 
##' @param features features in GRanges
##' @return list
##' @importFrom IRanges precede
##' @importFrom IRanges follow
##' @importFrom IRanges start
##' @importFrom IRanges end
##' @importFrom BiocGenerics strand
## @importMethodsFrom GenomicRanges strand
##' @author G Yu
getNearestFeatureIndicesAndDistances <- function(peaks, features) {
    ## peaks only conatin all peak records, in GRanges object
    ## feature is the annotation in GRanges object

    ## only keep start position based on strand
    start(features) <- end(features) <- ifelse(strand(features) == "+", start(features), end(features))

    ## add dummy NA feature for peaks that are at the last or first feature 
    ## suggested by Michael Kluge
    features.bak <- features
    seqlevels(features) <- c(seqlevels(features), "chrNA")
    dummy <- GRanges("chrNA", IRanges(1,1))
    dummy$tx_id <- -1
    dummy$tx_name <- "NA"
    features <- append(features, dummy)
    dummyID <- length(features)

    
    ## nearest from peak start
    ps.idx <- follow(peaks, features)
    
    ## nearest from peak end
    pe.idx <- precede(peaks, features)
    
    na.idx <- is.na(ps.idx) & is.na(pe.idx)
    ## if (sum(na.idx) > 1) {
    if (sum(na.idx) > 0) { ## suggested by Thomas Schwarzl
        ps.idx <- ps.idx[!na.idx]
        pe.idx <- pe.idx[!na.idx]
        peaks <- peaks[!na.idx]
    }

    
    # set NA values to dummy value if only one entry is affected
    ps.idx[is.na(ps.idx)] <- dummyID
    pe.idx[is.na(pe.idx)] <- dummyID

    
    ## features from nearest peak start
    psF <- features[ps.idx]
    ## feature distances from peak start
    ## psD <- ifelse(strand(psF) == "+",
    ##              start(peaks) - start(psF),
    ##              end(psF)-end(peaks))

    psD <- ifelse(strand(psF) == "+", 1, -1) *
        (start(peaks) - start(psF))
    
    ## features from nearest peak end
    peF <- features[pe.idx]
    ## feature distances from peak end
    ## peD <- ifelse(strand(peF) == "+",
    ##               end(peaks) - start(peF),
    ##               end(peF)-start(peaks))

    peD <- ifelse(strand(peF) == "+", 1, -1) *
        (end(peaks) - start(peF))

    psD[ps.idx == dummyID] <- Inf # ensure that there is even no match if a seq with name "chrNA" exists
    peD[pe.idx == dummyID] <- Inf # ensure that there is even no match if a seq with name "chrNA" exists
   
    pse <- data.frame(ps=psD, pe=peD)
    j <- apply(pse, 1, function(i) which.min(abs(i)))

    ## index
    idx <- ps.idx
    idx[j==2] <- pe.idx[j==2]
    
    ## distance
    dd <- psD
    dd[j==2] <- peD[j==2]

    
    hit <- findOverlaps(peaks, features)
    if ( length(hit) != 0 ) {
        qh <- queryHits(hit)
        hit.idx <- getFirstHitIndex(qh)
        hit <- hit[hit.idx]
        peakIdx <- queryHits(hit)
        featureIdx <- subjectHits(hit)

        idx[peakIdx] <- featureIdx
        dd[peakIdx] <- 0
    }

    ## pn.idx <- nearest(peaks, features)
    ## isOverlap <- sapply(1:length(pn.idx), function(i) {
    ##     isPeakFeatureOverlap(peaks[i], features2[pn.idx[i]])
    ## })
    ## isOverlap <- unlist(isOverlap)

    ## if(sum(isOverlap) > 0) {
    ##     idx[isOverlap] <- pn.idx[isOverlap]
    ##     dd[isOverlap] <- 0
    ## }
    
    res <- list(index=idx, distance=dd, peak=peaks)
    
    return(res)
}

isPeakFeatureOverlap <- function(peak, feature) {
    peakRange <- ranges(peak)
    featureRange <- ranges(feature)
    x <- intersect(peakRange, featureRange)
    return(length(x) != 0)
}
