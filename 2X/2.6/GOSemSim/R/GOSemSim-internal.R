ygcGetParents <- function(ont="MF") {
	if(!exists("GOSemSimEnv")) .initial()
	wh_Parents <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"	
	)		
	Parents <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFPARENTS) ,
		BP = AnnotationDbi::as.list(GOBPPARENTS) , 
		CC = AnnotationDbi::as.list(GOCCPARENTS)	
	)
	assign(eval(wh_Parents), Parents, envir=GOSemSimEnv)
}

ygcGetOffsprings <- function(ont="MF") {
	if(!exists("GOSemSimEnv")) .initial()
	wh_Offsprings <- switch(ont,
		MF = "MFOffsprings",
		BP = "BPOffsprings",
		CC = "CCOffsprings"	
	)		
	Offsprings <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFOFFSPRING) ,
		BP = AnnotationDbi::as.list(GOBPOFFSPRING) , 
		CC = AnnotationDbi::as.list(GOCCOFFSPRING)	
	)
	assign(eval(wh_Offsprings), Offsprings, envir=GOSemSimEnv)
}

ygcGetAncestors <- function(ont="MF") {
	if(!exists("GOSemSimEnv")) .initial()
	wh_Ancestors <- switch(ont,
		MF = "MFAncestors",
		BP = "BPAncestors",
		CC = "CCAncestors"	
	)		
	Ancestors <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFANCESTOR) ,
		BP = AnnotationDbi::as.list(GOBPANCESTOR) , 
		CC = AnnotationDbi::as.list(GOCCANCESTOR)	
	)
	assign(eval(wh_Ancestors), Ancestors, envir=GOSemSimEnv)
}

ygcCheckAnnotationPackage <- function(species){
	pkgname <- switch (species,
		human = "org.Hs.eg.db",
		fly = "org.Dm.eg.db",
		mouse = "org.Mm.eg.db",
		rat = "org.Rn.eg.db",
		yeast = "org.Sc.sgd.db",
		zebrafish = "org.Dr.eg.db",
		worm = "org.Ce.eg.db",
		arabidopsis = "org.At.tair.db",
		ecolik12 = "org.EcK12.eg.db", 
		bovine	= "org.Bt.eg.db",
		canine	= "org.Cf.eg.db", 
		anopheles	=	"org.Ag.eg.db", 
		ecsakai	=	"org.EcSakai.eg.db", 
		chicken	=	"org.Gg.eg.db", 
		chimp	=	"org.Pt.eg.db", 
		malaria	=	"org.Pf.plasmo.db", 
		rhesus	=	"org.Mmu.eg.db", 
		pig	= 	"org.Ss.eg.db", 
		xenopus	=	"org.Xl.eg.db"
	)
	p <- installed.packages()
	pn <- p[,1]
	if (sum(pn==pkgname) == 0) {
		print("The corresponding annotation package did not installed in this machine.")
		print("GOSemSim will install and load it automatically.")
		#source("http://bioconductor.org/biocLite.R")
		#biocLite(pkgname)
		install.packages(pkgname,repos="http://www.bioconductor.org/packages/2.6/data/annotation",type="source")
	}
	switch (species,
		human = library("org.Hs.eg.db"),
		fly = library("org.Dm.eg.db"),
		mouse = library("org.Mm.eg.db"),
		rat = library("org.Rn.eg.db"),
		yeast = library("org.Sc.sgd.db"),
		zebrafish = library("org.Dr.eg.db"),
		worm = library("org.Ce.eg.db"),
		arabidopsis = library("org.At.tair.db"),
		ecolik12 = library("org.EcK12.eg.db"),
		bovine	= library("org.Bt.eg.db"),
		canine	= library("org.Cf.eg.db"), 
		anopheles	=	library("org.Ag.eg.db"), 
		ecsakai	=	library("org.EcSakai.eg.db"), 
		chicken	=	library("org.Gg.eg.db"), 
		chimp	=	library("org.Pt.eg.db"), 
		malaria	=	library("org.Pf.plasmo.db"), 
		rhesus	=	library("org.Mmu.eg.db"), 
		pig	= library("org.Ss.eg.db"), 
		xenopus	=	library("org.Xl.eg.db")			
	)
}

ygcGetGOMap <- function(organism="human") {
	if(!exists("GOSemSimEnv")) .initial()
	ygcCheckAnnotationPackage(organism)
	species <- switch(organism,
		human = "Hs",
		fly = "Dm",
		mouse = "Mm",
		rat = "Rn",
		yeast = "Sc",
		zebrafish = "Dr",
		worm = "Ce",
		arabidopsis = "At",
		ecolik12 = "EcK12",
		bovine	= "Bt",
		canine	= "Cf", 
		anopheles	=	"Ag", 
		ecsakai	=	"EcSakai", 
		chicken	=	"Gg", 
		chimp	=	"Pt", 
		malaria	=	"Pf", 
		rhesus	=	"Mmu", 
		pig	= "Ss", 
		xenopus	=	"Xl"
	)
	gomap <- switch(organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecolik12 = org.EcK12.egGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		anopheles	=	org.Ag.egGO, 
		ecsakai	=	org.EcSakai.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		malaria	=	org.Pf.plasmoGO, 
		rhesus	=	org.Mmu.egGO, 
		pig	= org.Ss.egGO, 
		xenopus	=	org.Xl.egGO		
	)
	assign(eval(species), gomap, envir=GOSemSimEnv)
}

`ygcGetOnt` <-  function(gene, organism, ontology, dropCodes) {
	if(!exists("GOSemSimEnv")) .initial()
	species <- switch(organism,
		human = "Hs",
		fly = "Dm",
		mouse = "Mm",
		rat = "Rn",
		yeast = "Sc",
		zebrafish = "Dr",
		worm = "Ce",
		arabidopsis = "At",
		ecolik12 = "EcK12",
		bovine	= "Bt",
		canine	= "Cf", 
		anopheles	=	"Ag", 
		ecsakai	=	"EcSakai", 
		chicken	=	"Gg", 
		chimp	=	"Pt", 
		malaria	=	"Pf", 
		rhesus	=	"Mmu", 
		pig	= "Ss", 
		xenopus	=	"Xl"
	)
	if (!exists(species, envir=GOSemSimEnv)) {
		ygcGetGOMap(organism)
	}	
	gomap <- get(species, envir=GOSemSimEnv)
	
    allGO <- gomap[[gene]]
    if (is.null(allGO)) {
    	return (NA)
    }
    if (sum(!is.na(allGO)) == 0) {
    	return (NA)
    }
    if(!is.null(dropCodes)) { 
      evidence<-sapply(allGO, function(x) x$Evidence) 
      drop<-evidence %in% dropCodes
      allGO<-allGO[!drop]
    }
    category<-sapply(allGO, function(x) x$Ontology)
    allGO<-allGO[category %in% ontology]

    if(length(allGO)==0) return (NA)
    return (unlist(unique(names(allGO))))
}


ygcWangMethod <- function(GOID1, GOID2, ont="MF", organism="human") {
	if(!exists("GOSemSimEnv")) .initial()
	weight.isa = 0.8
	weight.partof = 0.6

	if (GOID1 == GOID2)
		return (gosim=1)		

	Parents.name <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"
	)
	if (!exists(Parents.name, envir=GOSemSimEnv)) {
		ygcGetParents(ont)
	}
	Parents <- get(Parents.name, envir=GOSemSimEnv)
	
	sv.a <- 1
	sv.b <- 1
	sw <- 1
	names(sv.a) <- GOID1
	names(sv.b) <- GOID2 
	
	sv.a <- ygcSemVal(GOID1, ont, Parents, sv.a, sw, weight.isa, weight.partof)
	sv.b <- ygcSemVal(GOID2, ont, Parents, sv.b, sw, weight.isa, weight.partof)
	
	sv.a <- uniqsv(sv.a)
	sv.b <- uniqsv(sv.b)
	
	idx <- intersect(names(sv.a), names(sv.b))
	inter.sva <- unlist(sv.a[idx])
	inter.svb <- unlist(sv.b[idx])
	sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
	return(sim)
}



uniqsv <- function(sv) {
	sv <- unlist(sv)
	una <- unique(names(sv))
	sv <- unlist(sapply(una, function(x) {max(sv[names(sv)==x])}))
	return (sv)
}

ygcSemVal_internal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	p <- Parents[goid]
	p <- unlist(p[[1]])
	if (length(p) == 0) {
		#warning(goid, " may not belong to Ontology ", ont)
		return(0)
	}
	relations <- names(p)
	old.w <- w
	for (i in 1:length(p)) {
		if (relations[i] == "is_a") {
			w <- old.w * weight.isa
		} else {
			w <- old.w * weight.partof
		}
		names(w) <- p[i]
		sv <- c(sv,w)
		if (p[i] != "all") {
			sv <- ygcSemVal_internal(p[i], ont, Parents, sv, w, weight.isa, weight.partof)
		}
	}
	return (sv)
}

ygcSemVal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	if(!exists("GOSemSimCache")) return(ygcSemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof))
	goid.ont <- paste(goid, ont, sep=".")
	if (!exists(goid.ont, envir=GOSemSimCache)) {
	  	value <- ygcSemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof)
	  	assign(goid.ont, value, envir=GOSemSimCache)
		#cat("recompute ", goid, value, "\n")
	}
	else{
		#cat("cache ", goid, get(goid, envir=GOSemSimCache), "\n")
	}
	return(get(goid.ont, envir=GOSemSimCache))
}

`ygcInfoContentMethod` <- function(GOID1, GOID2, ont, measure, organism) {
	if(!exists("GOSemSimEnv")) .initial()
	fname <- paste("Info_Contents", ont, organism, sep="_")
	tryCatch(utils::data(list=fname, package="GOSemSim", envir=GOSemSimEnv))
	Info.contents <- get("IC", envir=GOSemSimEnv)

	rootCount <- max(Info.contents[Info.contents != Inf])
	Info.contents["all"] = 0
	p1 <- Info.contents[GOID1]/rootCount
	p2 <- Info.contents[GOID2]/rootCount    

	if (p1 == 0 || p2 == 0) return (NA)
	Ancestor.name <- switch(ont,
		MF = "MFAncestors",
		BP = "BPAncestors",
		CC = "CCAncestors"	
	)	
	if (!exists(Ancestor.name, envir=GOSemSimEnv)) {
		ygcGetAncestors(ont)
	}
	
	Ancestor <- get(Ancestor.name, envir=GOSemSimEnv)					
	ancestor1 <- unlist(Ancestor[GOID1])
	ancestor2 <- unlist(Ancestor[GOID2])
	if (GOID1 == GOID2) {
		commonAncestor <- GOID1
	} else if (GOID1 %in% ancestor2) {
		commonAncestor <- GOID1
	} else if (GOID2 %in% ancestor1) {
		commonAncestor <- GOID2
	} else { 
		commonAncestor <- intersect(ancestor1, ancestor2)
	}
	if (length(commonAncestor) == 0) return (NA)
	pms <- max(Info.contents[commonAncestor], na.rm=TRUE)/rootCount
	sim<-switch(measure,
   	    Resnik = pms,
   	    Lin = pms/(p1+p2),
   	    Jiang = 1 - min(1, -2*pms + p1 + p2), 
   	    Rel = 2*pms/(p1+p2)*(1-exp(-pms*rootCount))
	)   	
	return (sim)
}



ygcCompute_Information_Content <- function(dropCodes="NULL", ont, organism) {
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
	ygcCheckAnnotationPackage(wh_organism)
	gomap <- switch(wh_organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecolik12 = org.EcK12.egGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		anopheles	=	org.Ag.egGO, 
		ecsakai	=	org.EcSakai.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		malaria	=	org.Pf.plasmoGO, 
		rhesus	=	org.Mmu.egGO, 
		pig	= org.Ss.egGO, 
		xenopus	=	org.Xl.egGO		
	)
	mapped_genes <- mappedkeys(gomap)
	gomap = AnnotationDbi::as.list(gomap[mapped_genes])
	if (!is.null(dropCodes)){
		gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology %in% wh_ont)))
		gomap<-sapply(gomap, function(x) x[2,x[1,]=="FALSE"])
		gomap<-gomap[sapply(gomap,length) >0]		
	}else {
		gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology %in% wh_ont))
	}
	
	goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE) # all GO terms appearing in an annotation	
	goids <- toTable(GOTERM)
	# all go terms which belong to the corresponding category..
	goids <- unique(goids[goids[,"Ontology"] == wh_ont, "go_id"])  	
	gocount <- table(goterms)
	goname <- names(gocount) #goid of specific organism and selected category.
	## ensure goterms not appearing in the specific annotation have 0 frequency..
	go.diff <- setdiff(goids, goname)
	m <- double(length(go.diff)) 
	names(m) <- go.diff
	gocount <- as.vector(gocount)
	names(gocount) <- goname
	gocount <- c(gocount, m)

	Offsprings.name <- switch(wh_ont,
		MF = "MFOffsprings",
		BP = "BPOffsprings",
		CC = "CCOffsprings"	
	)	
	if (!exists(Offsprings.name, envir=GOSemSimEnv)) {
		ygcGetOffsprings(wh_ont)
	}
	Offsprings <- get(Offsprings.name, envir=GOSemSimEnv)	
	cnt <- sapply(goids,function(x){ c=gocount[unlist(Offsprings[x])]; gocount[x]+sum(c[!is.na(c)])})		
	names(cnt) <- goids	
	IC<- -log(cnt/sum(gocount))		
	save(IC, file=paste(paste("Info_Contents", wh_ont, organism, sep="_"), ".rda", sep=""))
}

rebuildICdata <- function(){
	ont <- c("MF","CC", "BP")
	species <- c("human", "rat", "mouse", "fly", "yeast", "zebrafish", "arabidopsis","worm", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus") 
	cat("------------------------------------\n")
	cat("calulating Information Content...\nSpecies:\t\tOntology\n")
	for (i in ont) {
		for (j in species) {
			cat(j)
			cat("\t\t\t")
			cat(i)
			cat("\n")
			ygcCompute_Information_Content(ont=i, organism=j)
		}
	}
	cat("------------------------------------\n")
	print("done...")
}