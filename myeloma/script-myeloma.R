library(multicore)

## #####################################
## file of the main user-defined parameters
source("README_userDefinedParameters_myeloma.R")

## ##########################
## load function files from sourcedir above
source(paste(sourcedir,"usefulFunctions.R",sep="/")) ## data prep functions, etc
source(paste(sourcedir,"modelFunctions.R",sep="/")) ## load model calculation functions

## #########################
## scripts that loads the data
if( arrayplatform=="agilent" ){
  source(paste(sourcedir,"loaddata_Agilent.R",sep="/"))
}else if( arrayplatform=="affymetrix" ){
  source(paste(sourcedir,"loaddata_Affy.R",sep="/"))
}else if( arrayplatform=="GEO" ){
  source(paste(sourcedir,"loaddata_GEO.R",sep="/"))
}else if( arrayplatform=="expressiontable" ){
  source(paste(sourcedir,"load_expressiontable.R",sep="/"))
}else if( arrayplatform=="simulation" ){
  source(paste(sourcedir,"simulateMiRmRNAdata.R",sep="/"))
}

if( exists("typeadder.script") && !is.null(typeadder.script) ){
  source(typeadder.script)
}

if( exists("typefollows.generator") && !is.null(typefollows.generator) ){
  source(typefollows.generator)
}


## #########################################

## this script initializes some parameters and data, and reduces the
## data set so that we use only those miRs and mRNAs that: (1) are
## "present" with significance in a certain number of sample types,
## (2) have a certain signal-to-noise ratio, and (3) could participate
## in a target pairing according to the loaded database data. The
## script also re-scales the data for each gene to a mean of 0 and
## standard deviation of 1, saving the mean log expression and
## standard deviation into *possmirs* for possible influence on the
## targeting interactions
source(paste(sourcedir,"modelDataInitializations.R",sep="/"))


## ######################################
## MODEL FIT

doclustering <- list(mir=TRUE,m=FALSE)

initialiters <- 5
dataiters <- 1
smalldataiters <- 3
interactioniters <- 10
bigiters <- 1 #25

info$nclus$mir <- 15

doside <- "mir"
doclusnums <- rep(15,5)
reps <- 1
mc.cores <- 3

allret <- mclapply(doclusnums,function(nc){
  info$nclus[[doside]] <- nc
  ret2 <- lapply( 1:reps,function(i){
    ret <- fitModel(xpr,pars,priors,info,possmirs,
                    initialiters=initialiters,
                    dataiters=dataiters,
                    smalldataiters=smalldataiters,
                    interactioniters=interactioniters,
                    bigiters=bigiters,
                    doclustering=doclustering)
    return(ret)
  })
  return(ret2)
},mc.cores=mc.cores)


## save(allret,file="allret.Rdata")

## load("allret.Rdata")


best.ind <- which.max( sapply(allret,function(x)
                              x[[1]]$llh[[1]]) )
## best.ind <- 1

ret <- allret[[best.ind]][[1]]
pars <- ret$pars
priors <- ret$priors
info <- ret$info
llh <- ret$llh



## ################################################
## make a table of the probability of targeting and then sort it

source(paste(sourcedir,"generateNiceResults.R",sep="/"))

tablelist <- generateNiceResults(xpr,llh,
                                 pars,priors,
                                 info,
                                 possmirs)

write.csv(tablelist$intz,file="resultsTable.csv")

as.matrix(unique(tablelist$intz$mRNA.main[1:100]))



## #############################################
##



## ###################################################
## GO and KEGG analysis on the DEG by the SAM procedure

## {
##   source("http://bioconductor.org/biocLite.R")
##   biocLite("GOstats")
##   biocLite("multtest")
##   biocLite("gcrma")
##   biocLite("RankProd")
##   biocLite("siggenes")
##   biocLite("biomaRt")
## }

## install.packages("EMA")
## require(EMA)

## ################
## hgu133a is the array Lionetti used
##
## source("http://bioconductor.org/biocLite.R")
## biocLite("hgu133a.db")
require(hgu133a.db) ## this is automatic

## we need Entrez IDs; here is a translator
probenames <- info$genes$m
names(probenames) <- probenames
entrezids <- sapply(probenames,function(pn)
                    xpr$m$genes[pn,"Entrez ID"])

## this maps entrezIDs to probe names used in the "universe"
uniprobe2entrez <- unlist(as.list(hgu133aENTREZID))
entrez2uniprobe <- names(uniprobe2entrez)
names(entrez2uniprobe) <- uniprobe2entrez

## ##################
## remove duplicated probe sets
tablelist.nodups <- lapply(tablelist,function(tab){
    ind <- duplicated(tab[,c("probe.m","probe.mir")])
    tab.nodups <- tab[!ind,]
})


## ########################
## the lists used for KEGG analysis
probes.forKEGG <- lapply( c(100,200),function(topn){
  lapply(tablelist.nodups,function(tab){
    ##lapply(tl,function(tab){
      test.entrezids <- entrezids[tab$probe.m[1:topn]]
      test.probes <- entrez2uniprobe[as.character(test.entrezids)]
      return(test.probes)
    ##})
  })
})



temp <- table(names(probes.forKEGG[[1]]$intz))
top100tab <- cbind(entrez2uniprobe[names(temp)],
                   names(temp),
                   sapply(names(temp),function(x)
                          xpr$m$genes$ORF[xpr$m$genes$`Entrez ID`==x]),
                   temp)[order(entrez2uniprobe[names(temp)]),]
rownames(top100tab) <- NULL
colnames(top100tab) <- c("probe","entrez","name","occurrences")


write.table(top100tab,file=("top100.csv"),
            sep="\t",quote=FALSE,row.names=FALSE)


## write.table(probes.forGOandKEGG[[2]]$IO$intz,
##             file="probesForGOandKEGG_200_IO_intz.txt",
##             col.names=FALSE,row.names=FALSE,quote=FALSE)
