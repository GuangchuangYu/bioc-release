##' measuring similarities bewteen two gene vectors.
##'
##' provide two entrez gene vectors, this function will calculate their similarity.
##' @title geneSim
##' @param geneID1 entrez gene vector
##' @param geneID2 entrez gene vector
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param combine One of "max", "average", "rcmax", "BMA" methods, for combining semantic similarity scores of multiple DO terms associated with gene/protein.
##' @return score matrix
##' @importFrom DO.db DOPARENTS
##' @importFrom DO.db DOANCESTOR
##' @importFrom GOSemSim combineScores
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
geneSim <- function(geneID1,
                    geneID2,
                    measure="Wang",
                    combine="BMA") {

    DOID1 <- lapply(geneID1, gene2DO)
    DOID2 <- lapply(geneID2, gene2DO)
    m <- length(geneID1)
    n <- length(geneID2)
    scores <- matrix(NA, nrow=m, ncol=n)
    rownames(scores) <- geneID1
    colnames(scores) <- geneID2

    for (i in 1:m) {
        for (j in 1:n) {
            if(any(!is.na(DOID1[[i]])) &&  any(!is.na(DOID2[[j]]))) {
                s <- doSim(DOID1[[i]],
                           DOID2[[j]],
                           measure
                           )
                scores[i,j] = combineScores(s, combine)
            }
        }
    }
    if (nrow(scores) == 1 & ncol(scores) == 1)
        scores = as.numeric(scores)
    return(scores)
}
