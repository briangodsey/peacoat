## UPDATED: miRanda predictions also added

## this script loads the TargetScan data file and reduces the table to
## only relevant genes (*ts.cons* for conserved sites) and a list of
## possible targeting miRs, listed by mRNA names (*possmirs.cons*)

## ################################################
## read TargetScan file (downloaded from the web site)
## and
## miRanda files (downloaded from web site)

if( data.organism=="human" ){
  ## HUMAN
  ts.cons <- read.delim(paste(databasefiledir,"predicted/targetscan/human/Conserved_Site_Context_Scores.txt",
                              sep="/"),stringsAsFactors=FALSE)
  miRandadat <- read.delim(paste(databasefiledir,"predicted/miRanda/hg19_predictions_S_C_aug2010.txt",
                           sep="/"),stringsAsFactors=FALSE)
}else if( data.organism=="mouse" ){
  ## MOUSE
  ts.cons <- read.delim(paste(databasefiledir,"predicted/targetscan/mouse/Conserved_Site_Context_Scores.txt",
                        sep="/"),stringsAsFactors=FALSE)
  miRandadat <- read.delim(paste(databasefiledir,"predicted/miRanda/mm9_predictions_S_C_aug2010.txt",
                           sep="/"),stringsAsFactors=FALSE)
}else{
  cat("There is no database data for your organism.\n")
  ts.cons <- data.frame(miRNA="0",Gene.Symbol="0",stringsAsFactors=FALSE)
  miRandadat <- data.frame(mirna_name="0",gene_symbol="0",stringsAsFactors=FALSE)
}

## ## gsymbs <- intersect(unique(ts.cons$Gene.Symbol),
## ##                     fit.m$genes$GeneName)
## gsymbs <- unique(ts.cons$Gene.Symbol)
## names(gsymbs) <- gsymbs


## ###################################################


## if expression data are already loaded, reduce tables to those genes
ind <- which(ts.cons[,"miRNA"] %in%
             unique(xpr.orig$mir$genes[[info$possmir.gname$mir]]))
ts.cons <- ts.cons[ind,]

ind <- which(ts.cons[,"Gene.Symbol"] %in%
             unique(xpr.orig$m$genes[[info$possmir.gname$m]]))
ts.cons <- ts.cons[ind,]

ind <- which(miRandadat[,"mirna_name"] %in%
             unique(xpr.orig$mir$genes[[info$possmir.gname$mir]]))
miRandadat <- miRandadat[ind,]

ind <- which(miRandadat[,"gene_symbol"] %in%
             unique(xpr.orig$m$genes[[info$possmir.gname$m]]))
miRandadat <- miRandadat[ind,]

## ###################################
## find and remove duplicate rows

## n <- dim(ts.cons)[1]
## ab <- ts.cons[1:(n-1),]
## ac <- ts.cons[2:n,]
## ind <- which(apply(ab==ac,1,all))

## there's a function for that, dummy
dups <- !duplicated(ts.cons[, c("Gene.Symbol","miRNA",
                                "UTR_start","UTR.end") ])
ts.cons <- ts.cons[dups,]

dups <- !duplicated(miRandadat[, c("mirna_name","gene_symbol",
                                   "gene_start","gene_end") ])
miRandadat <- miRandadat[dups,]



## ######################################
## group possible targeting miRs by target mRNA

## for ts.cons
## takes a couple of minutes
if( nrow(ts.cons)>0 ){
  gscol <- ts.cons$Gene.Symbol
  possmirs.tscons <- list()
  istart <- 1
  curgs <- gscol[1]
  for(i in 1:dim(ts.cons)[1]){ # 100000){
    if( curgs!=gscol[i] || i==dim(ts.cons)[1] ){
      possmirs.tscons[[curgs]] <- ts.cons[istart:(i-1),]
      curgs <- gscol[i]
      istart <- i
    }
  }
}else{
  possmirs.tscons <- list()
}


## for miRandadat
## takes a couple of minutes
if( nrow(miRandadat)>0 ){
  gscol <- miRandadat$gene_symbol
  possmirs.miRanda <- list()
  istart <- 1
  curgs <- gscol[1]
  for(i in 1:dim(miRandadat)[1]){ # 100000){
    if( curgs!=gscol[i] || i==dim(miRandadat)[1] ){
      possmirs.miRanda[[curgs]] <- miRandadat[istart:(i-1),]
      curgs <- gscol[i]
      istart <- i
    }
  }
}else{
  possmirs.miRanda <- list()
}


## #########################################
## combine the possmirs.* into one list
tsnames <- c(miRNA.main="miRNA",mRNA.main="Gene.Symbol")
mirandanames <- c(miRNA.main="mirna_name",mRNA.main="gene_symbol")

## allgenes <- union( names(possmirs.tscons),
##                   names(possmirs.miRanda) )
allgenes <- unique(xpr.orig$m$genes[,info$possmir.gname$m])
names(allgenes) <- allgenes

possmirs.orig <- lapply(allgenes,function(gn){
  ##cat(gn,"\n")
  tabnames <- c("miRNA.main","mRNA.main",
                "TS.pred","miRanda.pred",
                tsnames,colnames(ts.cons),
                mirandanames,colnames(miRandadat))
  names(tabnames) <- tabnames

  if( length(possmirs.tscons)>0 && gn %in% names(possmirs.tscons) ){
    tsmat <- possmirs.tscons[[gn]]
    if( is.null(tsmat) ){
      tsmat.aug <- as.data.frame(t(sapply(tabnames,function(x) 1 ))[0,])
      ##tsmat <- possmirs.tscons[[1]][0,,drop=FALSE]
    }else{
      tsmat.aug <- sapply(tabnames,function(tn){
        if( tn %in% c("miRNA.main","mRNA.main") ){
          ret <- tsmat[,tsnames[tn]]
        }else if( tn == "TS.pred"){
          ret <- rep(1,dim(tsmat)[1])
        }else if( tn %in% names(ts.cons) ){
          ret <- tsmat[,tn]
        }else{
          ret <- rep(0,dim(tsmat)[1])
        }
        return(ret)
      })
      if( is.null(dim(tsmat.aug)) ) tsmat.aug <- t(tsmat.aug)
    }      
  }else{
    tsmat.aug <- as.data.frame(t(sapply(tabnames,function(x) 1 ))[0,])
  }

  if( length(possmirs.miRanda)>0 ){
    mirandamat <- possmirs.miRanda[[gn]]
    if( is.null(mirandamat) ){
      mirandamat.aug <- as.data.frame(t(sapply(tabnames,function(x) 1 ))[0,])
      ##mirandamat <- possmirs.miRanda[[1]][0,,drop=FALSE]
    }else{
      mirandamat.aug <- sapply(tabnames,function(tn){
        if( tn %in% c("miRNA.main","mRNA.main") ){
          ret <- mirandamat[,mirandanames[tn]]
        }else if( tn == "miRanda.pred"){
          ret <- rep(1,dim(mirandamat)[1])
        }else if( tn %in% names(miRandadat) ){
          ret <- mirandamat[,tn]
        }else{
          ret <- rep(0,dim(mirandamat)[1])
        }
        return(ret)
      })
      if( is.null(dim(mirandamat.aug)) ) mirandamat.aug <- t(mirandamat.aug)
    }
  }else{
    mirandamat.aug <- as.data.frame(t(sapply(tabnames,function(x) 1 ))[0,])
  }

  if( dim(tsmat.aug)[1]==0 ){
    retmat <- mirandamat.aug
  }else if( dim(mirandamat.aug)[1]==0 ){
    retmat <- tsmat.aug
  }else{
    retmat <- rbind(tsmat.aug,
                    mirandamat.aug)
  }
  rownames(retmat) <- NULL
  
  return(retmat)
})
                      
if( data.organism != "simulation" ){
  save(possmirs.orig,
       file=paste(generatedfiledir,"possmirs_orig.Rdata",sep="/"))
}

## save(ts.cons,file="tscons.Rdata")
## save(miRandadat,file="miRandadat.Rdata")


