setClass("candidates", representation(bamFilePath = "character",
                                      candidatesFileName= "character",
                                      build="character",
                                      readLength = "numeric",
                                      mmRate = "numeric",
                                      indelRate = "numeric",
                                      quickScore = "numeric",
                                      fullScore = "numeric",
                                      forplot = "list"))

setGeneric("bamFilePath", function(object) standardGeneric("bamFilePath"))
setGeneric("bamFilePath<-", function(object,value) standardGeneric("bamFilePath<-"))
setReplaceMethod("bamFilePath","candidates", function(object,value){object@bamFilePath<- value;object})


setGeneric("candidatesFileName", function(object) standardGeneric("candidatesFileName"))
setGeneric("candidatesFileName<-", function(object,value) standardGeneric("candidatesFileName<-"))
setReplaceMethod("candidatesFileName","candidates", function(object,value){object@candidatesFileName<- value;object})



setGeneric("build", function(object) standardGeneric("build"))
setGeneric("build<-", function(object,value) standardGeneric("build<-"))
setReplaceMethod("build","candidates", function(object,value){object@build<- value;object})


setGeneric("readLength", function(object) standardGeneric("readLength"))
setGeneric("readLength<-", function(object,value) standardGeneric("readLength<-"))
setReplaceMethod("readLength","candidates", function(object,value){object@readLength<- value;object})


setGeneric("mmRate", function(object) standardGeneric("mmRate"))
setGeneric("mmRate<-", function(object,value) standardGeneric("mmRate<-"))
setReplaceMethod("mmRate","candidates", function(object,value){object@mmRate<- value;object})


setGeneric("indelRate", function(object) standardGeneric("indelRate"))
setGeneric("indelRate<-", function(object,value) standardGeneric("indelRate<-"))
setReplaceMethod("indelRate","candidates", function(object,value){object@indelRate<- value;object})



setGeneric("quickScore", function(object,...) standardGeneric("quickScore"))
#setGeneric("quickScore<-", function(object,...) standardGeneric("quickScore<-"))
setMethod("quickScore","candidates", function(object,...){

    object@quickScore<-doQuickScore(filename=object@candidatesFileName,
                                  readLength=object@readLength,
                                  bamFilePath=object@bamFilePath,
                                  mmRate=object@mmRate,
                                  indelRate=object@indelRate,
                                  ...)
    object
})



 setGeneric("fullScoreAndView", function(object,...) standardGeneric("fullScoreAndView"))
# setGeneric("fullScoreAndView<-", function(object,value) standardGeneric("fullScoreAndView<-"))
 setMethod("fullScoreAndView","candidates", function(object,...){
     retvals <-ViewAndScoreFull(filename=object@candidatesFileName,
                                  readLength=object@readLength,
                                  bamFilePath=object@bamFilePath,
                                  mmRate=object@mmRate,
                                  indelRate=object@indelRate,
                                  build=object@build,
                                ...)

     object@forplot <- retvals[['forplot']]
     object@fullScore <- retvals[['score']]
     object
 })



 setGeneric("plotSV", function(object,indices=0,junctionsToShow=c('sv','side1','side2','all')[4],pdfname='',
                               width=8,height=7,config='TwoRows',...) standardGeneric("plotSV"))
 setMethod("plotSV","candidates", function(object,indices=0,junctionsToShow=c('sv','side1','side2','all')[4],pdfname='',width=8,height=7,
                                           config='TwoRows',...){
     if(indices[1]==0) indices <- seq_along(object@forplot)

     if(pdfname != '') pdf(pdfname,width=width,height=height)
     for(ii in indices){
     svPlot <-formatPlot(dfplot=object@forplot[[ii]][[1]],...)
     if(junctionsToShow != 'sv'){
     contig1Plot <- formatPlot(dfplot=object@forplot[[ii]][[2]],plotJunctions=svPlot[[2]],nolegend=TRUE,...)
     contig2Plot <- formatPlot(dfplot=object@forplot[[ii]][[3]],plotJunctions=svPlot[[2]],nolegend=TRUE,...)
     columnprintContig1<-which.min(abs(as.numeric(svPlot[[2]])-as.numeric(c(contig1Plot[[2]],contig2Plot[[2]]))))
     if(columnprintContig1==2){
         temp1<-contig1Plot
         contig1Plot <- contig2Plot
         contig2Plot <- temp1
     }

 }

     if(pdfname=='') dev.new()

     if(junctionsToShow=='sv') print(svPlot[[1]])
     else if (junctionsToShow=='side1') print(contig1Plot[[1]])
     else if (junctionsToShow=='side2') print (contig2Plot[[1]])
     else {
         if(config=='TwoRows'){
         grid.newpage()
         pushViewport(viewport(layout = grid.layout(2, 1,heights = unit(c(.55,.45),"null"))))
         print(svPlot[[1]],vp = viewport(layout.pos.row = 1, layout.pos.col=1))
         pushViewport(viewport(layout = grid.layout(2, 2,heights=unit(c(.55,.45),"null"))))
         print(contig1Plot[[1]],vp = viewport(layout.pos.row = 2, layout.pos.col=1))
         print(contig2Plot[[1]],vp = viewport(layout.pos.row = 2, layout.pos.col=2))
     }else if(config=='ThreeRows'){
         grid.newpage()
         pushViewport(viewport(layout = grid.layout(3, 1)))
         print(svPlot[[1]],vp = viewport(layout.pos.row = 1, layout.pos.col=1))
         print(contig1Plot[[1]],vp = viewport(layout.pos.row = 2, layout.pos.col=1))
         print(contig2Plot[[1]],vp = viewport(layout.pos.row = 3, layout.pos.col=1))
     }

     }

     cat('.')
 }
     if(pdfname != '') graphics.off()
 })

## setGeneric("ViewSV", function(object,...) standardGeneric("ViewSV"))
## setMethod("ViewSV",signature(object="candidates"),
##           function(object,indices,...){
##               ## do plotting
##           })



## setMethod("gcCorrect", signature(object="matrix"), function(object, ...){
##     gcCorrectMain(object, ...)
## })

## gcCorrectBeadStudioSet <- function(object, ...){
## 	args <- list(...)
## 	if("returnOnlyTV" %in% names(args)){
## 		is.score <- args[["returnOnlyTV"]]
## 	} else is.score <- FALSE
## 	r <- lrr(object)
## 	if(is(r, "matrix")){
## 		r <- r/100
## 	}
## 	isff <- is(r, "ff_matrix") ## need to be careful
## 	if(isff){
## 		index.list <- split(seq_len(ncol(object)), ocSamples())
## 	} else index.list <- list(seq_len(ncol(object)))
## 	if(is.score) score.list <- list()
## 	for(i in seq_along(index.list)){
## 		j <- index.list[[i]]
## 		pos <- position(object)
## 		chr <- paste("chr", chromosome(object), sep="")
## 		res <- gcCorrectMain(Ms=r,
## 				     chr=chr,
## 				     starts=pos,
## 				     samplechr=unique(chr),
## 				     build=genomeBuild(object),
## 				     ...)
## 		if(!is.score){  ## if not returning TV score, update the brList object
## 			res <- integerMatrix(res, 100)
## 			lrr(object) <- res
## 		} else {
## 			score.list[[i]] <- res
## 		}
## 	}
## 	if(is.score) {
## 		results <- score.list
## 	} else results <- object
## 	return(results)
## }

## setMethod("gcCorrect", signature(object="BeadStudioSet"),
## 	  function(object, ...){
## 		  gcCorrectBeadStudioSet(object, ...)
## 	  })
## setMethod("gcCorrect", signature(object="BafLrrSet"),
## 	  function(object, ...){
## 		  gcCorrectBeadStudioSet(object, ...)
## 	  })



## gcCorrectBafLrrList <- function(object, index.samples, providedGC=NULL,...){
## 	args <- list(...)
## 	if("returnOnlyTV" %in% names(args)){
## 		return.score <- args[["returnOnlyTV"]]
## 	} else return.score <- FALSE
## 	if(return.score) stop(paste("return.score not implemented for ", class(object), " objects"))
## 	r <- lrr(object)
## 	isff <- is(r[[1]], "ff")
## 	if(missing(index.samples))
## 		index.samples <- seq_len(ncol(object[[1]]))
## 	## to keep RAM in check, do in batches of samples
## 	index.list <- splitIndicesByLength(index.samples, ocSamples())
## 	l <- elementLengths(object)
## 	chr <- paste("chr", rep(chromosome(object), l), sep="")
## 	pos <- unlist(position(object))
## 	##if(return.score) score.list <- list()
## 	.packages <- c("oligoClasses", "ArrayTV")
## 	isFFloaded <- isPackageLoaded("ff")
## 	if(isFFloaded) .packages <- c("ff", .packages)
## 	gc.provided <- !is.null(providedGC)
## 	if(gc.provided) providedGC <- as.matrix(providedGC)
## 	if("returnOnlyTV" %in% names(list(...))){
## 		only.tv <- list(...)[["returnOnlyTV"]]
## 	} else only.tv <- FALSE
## 	j <- 1
## 	reslist <- foreach(j=index.list, .packages=.packages) %dopar% {
## 		rr <- lapply(r, function(x, j) x[, j, drop=FALSE]/100, j=j)
## 		R <- do.call(rbind, rr)
## 		rm(rr)
## 		namatrix <- is.na(R)
## 		R[namatrix] <- 0 ## not ideal
## 		if(gc.provided){
## 			for(k in seq_len(ncol(providedGC))){
## 				R <- (gcCorrectMain(Ms=R,
## 						   chr=chr,
## 						   starts=pos,
## 						   samplechr=unique(chr),
## 						   build=genomeBuild(object),
## 						   providedGC=providedGC[, k],
## 						   ...))[['correctedVals']]
## 			}
## 		} else {
## 			R <- (gcCorrectMain(Ms=R,
## 					   chr=chr,
## 					   starts=pos,
## 					   samplechr=unique(chr),
## 					   build=genomeBuild(object),
## 					   providedGC=providedGC,
## 					   ...))[['correctedVals']]
## 		}
## 		if(!only.tv) {
## 			R <- integerMatrix(R, 100)
## 			R[namatrix] <- NA
## 		}
## 		##rm(R); gc()
## 		##
## 		## update ff object.  Writing is expensive -- only do
## 		## once
## 		##
## 		## Only reasonable to do this when ff package is
## 		## loaded. Otherwise, the replacment method is
## 		## transient and we do not want to return an entire
## 		## copy of the object
## 		if(isFFloaded & !only.tv) {
## 			lrr(object) <- R
## 			R <- NULL
## 		}
## 		gc()
## 		return(R)
## 	}
## 	if(!only.tv){
## 		if(!isFFloaded){
## 			res <- do.call("cbind", reslist)
## 			lrr(object) <- res
## 		} else {
## 			## we do not want to return an entire copy of the
## 			## object if ff package is loaded
## 			object <- NULL
## 		}
## 	} else object <- reslist
## 	return(object)
## }

## setMethod("gcCorrect", signature(object="BafLrrSetList"),
## 	  function(object, ...){
## 		  gcCorrectBafLrrList(object, ...)
## 	  })

## gcModel <- function(data, window, verbose=FALSE){
## 	## assume data is summarized experiment
## 	build <- metadata(rowData(data))$genome
## 	library(paste("BSgenome.Hsapiens.UCSC.", build, sep=''), character.only=TRUE)
## 	if(verbose) print('Getting gc content From BS genome Object')
## 	Hsapiens <- get("Hsapiens")
## 	chroms <- unique(chromosome(data))
## 	maxwin <- increm <- window
## 	## query genome prior to looking at genomic position of
## 	## markers??
## 	gc <- list()
## 	for(i in seq_along(chroms)){
## 		if(verbose) cat(".")
## 		chr <- chroms[i]
## 		j <- which(chromosome(data)==chr)
## 		## calculates the gc content for each window of size increm
## 		pregcFrac <- letterFrequencyInSlidingView(unmasked(Hsapiens[[chr]]), view.width=increm, 'CG', as.prob=TRUE)
## 		##if(verbose) print('gc content stored')
## 		startinds <- rep(start(data)[j], each=maxwin/increm) +
## 			rep(seq(0, maxwin-increm, increm), length(j))
## 		startinds[which(startinds<1)] <- 1
## 		startinds[which(startinds>length(pregcFrac))] <- length(pregcFrac)
## 		startindbackwards <- rep(start(data)[j], each=maxwin/increm)-rep(seq(increm, maxwin, increm), length(j))
## 		gcFrac1 <- pregcFrac[startinds]
## 		if(any(is.na(gcFrac1)))
## 			gcFrac1[is.na(gcFrac1)] <- mean(gcFrac1, na.rm=TRUE)
## 		startindbackwards[which(startindbackwards<1)] <- 1
## 		startindbackwards[which(startindbackwards>length(pregcFrac))] <- length(pregcFrac)
## 		gcFracbackwards1 <- pregcFrac[startindbackwards]
## 		if(any(is.na(gcFracbackwards1))){
## 			gcFracbackwards1[is.na(gcFracbackwards1)] <- mean(gcFracbackwards1, na.rm=TRUE)
## 		}
## 		gc[[i]] <- (gcFrac1+gcFracbackwards1)/2
## 	}
## 	result <- unlist(gc)
## 	return(result)
## }

## midpoint <- function(object) start(object) + floor((width(object)-1)/2)

## .rescaleGC <- function(gc, cn, shift=0){
## 	x <- scale(gc) ## mean 0, sd1
## 	## give it same sd as cn
## 	x <- x*mad(cn,na.rm=TRUE)
## 	## give it same location as cn
## 	x <- x+mean(cn,na.rm=TRUE) + shift
## 	x
## }

## gcModelSeq <- function(data, window, verbose=FALSE){
## 	## assume data is summarized experiment
## 	build <- genome(rowData(data))[[1]]
## 	library(paste("BSgenome.Hsapiens.UCSC.", build, sep=''), character.only=TRUE)
## 	if(verbose) print('Getting gc content From BS genome Object')
## 	Hsapiens <- get("Hsapiens")
## 	chroms <- unique(chromosome(data))
## 	maxwin <- increm <- window
## 	starts <- midpoint(data)
## 	## query genome prior to looking at genomic position of
## 	## markers??
## 	gc <- list()
## 	for(i in seq_along(chroms)){
## 		if(verbose) cat(".")
## 		chr <- chroms[i]
## 		j <- which(chromosome(data)==chr)
## 		## calculates the gc content for each window of size increm
## 		pregcFrac <- letterFrequencyInSlidingView(unmasked(Hsapiens[[chr]]), view.width=increm, 'CG', as.prob=TRUE)
## 		##if(verbose) print('gc content stored')
## 		startinds <- rep(starts, each=maxwin/increm) +
## 			rep(seq(0, maxwin-increm, increm), length(j))
## 		startinds[which(startinds<1)] <- 1
## 		startinds[which(startinds>length(pregcFrac))] <- length(pregcFrac)
## 		startindbackwards <- rep(starts, each=maxwin/increm)-rep(seq(increm, maxwin, increm), length(j))
## 		gcFrac1 <- pregcFrac[startinds]
## 		if(any(is.na(gcFrac1)))
## 			gcFrac1[is.na(gcFrac1)] <- mean(gcFrac1, na.rm=TRUE)
## 		startindbackwards[which(startindbackwards<1)] <- 1
## 		startindbackwards[which(startindbackwards>length(pregcFrac))] <- length(pregcFrac)
## 		gcFracbackwards1 <- pregcFrac[startindbackwards]
## 		if(any(is.na(gcFracbackwards1))){
## 			gcFracbackwards1[is.na(gcFracbackwards1)] <- mean(gcFracbackwards1, na.rm=TRUE)
## 		}
## 		gc[[i]] <- (gcFrac1+gcFracbackwards1)/2
## 	}
## 	result <- unlist(gc)
## 	return(result)
## }

