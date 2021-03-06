##'Semantic Similarity Between two Genes
##'
##'Given two genes, this function will calculate the semantic similarity between
##'them, and return their semantic similarity and the corresponding GO terms
##'
##'
##'@param gene1 Entrez gene id.
##'@param gene2 Another entrez gene id.
##'@param ont One of "MF", "BP", and "CC" subontologies.
##'@param organism One of "anopheles", "arabidopsis", "bovine", "canine",
##'"chicken", "chimp", "coelicolor", "ecolik12","ecsakai", "fly", "human",
##'"malaria", "mouse", "pig", "rat","rhesus", "worm", "xenopus", "yeast" and
##'"zebrafish".
##'@param measure One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.
##'@param drop A set of evidence codes based on which certain annotations are
##'dropped. Use NULL to keep all GO annotations.
##'@param combine One of "max", "average", "rcmax", "BMA" methods, for combining
##'semantic similarity scores of multiple GO terms associated with protein or
##'multiple proteins assiciated with protein cluster.
##'@return list of similarity value and corresponding GO.
##'@seealso \code{\link{goSim}} \code{\link{mgoSim}} \code{\link{mgeneSim}}
##'\code{\link{clusterSim}} \code{\link{mclusterSim}}
##'@references Yu et al. (2010) GOSemSim: an R package for measuring semantic
##'similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
##'England), 26:7 976--978, April 2010. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
##'PMID: 20179076
##'
##'@keywords manip
##' @export
##'@examples
##'
##'	geneSim("241", "251", ont="MF", organism="human", measure="Wang")
##'
geneSim <- function(gene1, gene2, ont="MF", organism="human", measure="Wang", drop="IEA", combine="BMA"){
    go1 <- gene2GO(gene1, organism=organism, ont=ont, dropCodes=drop)
    go2 <- gene2GO(gene2, organism=organism, ont=ont, dropCodes=drop)
    if (is.na(go1) || is.na(go2))
        return (NA)
    res <- mgoSim(go1, go2, ont=ont, organism=organism, measure=measure, combine=combine)
    return (list(geneSim=res, GO1=go1, GO2=go2))
}

