`clusterSim` <- 
function(cluster1, cluster2, ont="MF", organism="human", measure="Wang", drop="IEA"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	
	size1 <- length(cluster1)
	size2 <- length(cluster2)
	if (size1 == 0 || size2 == 0) {
		return (NA)
	}
	gos1 <- lapply(cluster1, function(x) ygcGetOnt(x, organism= wh_organism, ontology= wh_ont, dropCodes=drop))
	gos2 <- lapply(cluster2, function(x) ygcGetOnt(x, organism= wh_organism, ontology= wh_ont, dropCodes=drop))

	allSim <- matrix(data=NA, nrow=size1, ncol=size2)
	assign("GOSemSimCache", new.env(hash=TRUE),envir=.GlobalEnv)
	for (i in 1:size1) {
		for (j in 1:size2){
			#allSim[i,j] <- geneSim(cluster1[i], cluster2[j], wh_ont, wh_organism, wh_measure, drop)$geneSim

			allSim[i,j] <- NA
			if(any(!is.na(gos1[i])) &&  any(!is.na(gos2[j])))
			{
				sim <- mgoSim(gos1[i],gos2[j], wh_ont, wh_organism, wh_measure)
				sim <- round(sim, digits=3)
				allSim[i,j] <- sim
			}
		}
	}	
	
	remove("GOSemSimCache", envir=.GlobalEnv)
	if (!sum(!is.na(allSim))) {
		return (NA)
	}
	
	result <- sum(allSim, na.rm=TRUE)/sum(!is.na(allSim))
	return (round(result, digits=3))

}