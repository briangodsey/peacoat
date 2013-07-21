library(limma)
library(affy)
library(gcrma)

## source("http://bioconductor.org/biocLite.R")
## biocLite("annotate")
## biocLite("hugene10stprobeset.db")
## biocLite("hugene10sttranscriptcluster.db")
## biocLite("mirna10probe")

library(annotate)
##library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)
##library(mirna10probe)


## TableOfSamplesFilename <- "~/mirtargeting/srujana_filenames_ALL.txt"

sampletable <- read.table(TableOfSamplesFilename,
                          sep="\t",header=TRUE,stringsAsFactor=FALSE)

fullfns <- list()
fullfns$m <- sampletable$filename[sampletable$RNAtype=="m"]
fullfns$mir <- sampletable$filename[sampletable$RNAtype=="mir"]

sampletypes <- list()
sampletypes$m <- sampletable$type[sampletable$RNAtype=="m"]
sampletypes$mir <- sampletable$type[sampletable$RNAtype=="mir"]


dat <- list()
dat.exprs <- list()
xpr.orig <- list()
for( side in c("mir","m") ){
  dat[[side]] <- lapply(fullfns[[side]],function(x){
    datx <- ReadAffy(filenames=x,cdfname=NULL)
    vals <- rma(datx,normalize=FALSE)
    return(vals)
  })

  dat.exprs[[side]] <- sapply(dat[[side]],exprs)

  probeids <- featureNames(dat[[side]][[1]])

  
  rownames(dat.exprs[[side]]) <- probeids

  xpr.orig[[side]] <- list()
  xpr.orig[[side]]$E <- normalizeBetweenArrays(dat.exprs[[side]],
                                               method="quantile")
  colnames(xpr.orig[[side]]$E) <- sampletypes[[side]]

  xpr.orig[[side]]$targets <- sampletable[sampletable$RNAtype==side,]

  if( side=="mir" ){
    tmp <- gsub("_st","",probeids) 
    probenames <- gsub("-star","*",tmp,fixed=TRUE)
    xpr.orig[[side]]$genes <- data.frame(ProbeName=probeids,
                                         GeneName=probenames,
                                         stringsAsFactors=FALSE)
    if( data.organism=="human" ){
      ind <- grep("hsa-",probenames,fixed=TRUE)
    }else if( data.organism=="mouse" ){
      ind <- grep("mmu-",probenames,fixed=TRUE)
    }
    xpr.orig[[side]] <- new("EList",xpr.orig[[side]])
    xpr.orig[[side]] <- xpr.orig[[side]][ind,]     
  }else if( side=="m" ){
    genenames <- getSYMBOL(probeids,"hugene10sttranscriptcluster.db")
    xpr.orig[[side]]$genes <- data.frame(ProbeName=probeids,
                                         GeneName=genenames,
                                         stringsAsFactors=FALSE)
    xpr.orig[[side]] <- new("EList",xpr.orig[[side]])
    xpr.orig[[side]] <- xpr.orig[[side]][!is.na(genenames),]     
  }
  

  gname[[side]] <- colnames(xpr.orig[[side]]$genes)[1]
  possmir.gname[[side]] <- colnames(xpr.orig[[side]]$genes)[2]
}




## for( side in c("mir","m"){
##   write.csv(dat.exprs.qu[[side]],
##             file=paste("affy_exprmat_",side,".txt",sep=""),
##             quote=FALSE,row.names=FALSE)
## }




