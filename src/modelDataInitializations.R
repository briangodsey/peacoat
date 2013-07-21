## #########################################
## various variable intitializations


## ## old parameters; changing not recommended
## onlygeneswithvalidation <- FALSE
## SNRthresh <- 0


## #######################################
## define some parameters
info <- NULL
info$nclus$mir <- 2 # temporary values
info$nclus$m <- 2

## column name in xpr[[side]]$genes that we want to use as main list
## of genes
##info$gname <- list(mir="SystematicName",m="GeneName")
info$gname <- gname

## column name in xpr[[side]]$genes that contains the names that
## *possmir* (database) uses
info$possmir.gname <- possmir.gname


## use these factors in determining interaction potential
info$interactionfactors <- c("constant",
                             "TS.pred","miRanda.pred",
                             ##"meanlogexp.m","meanlogexp.mir",
                             ##"maxlogexp.mir",
                             ##"Site.Type","UTR_start","UTR.end","X3prime.pairing",
                             ##"local.AU","position","TA","SPS",
                             ##"mirGCcontent",
                             ##"miRlength",
                             "context..score",
                             ##"context..score.percentile",
                             "mirsvr_score"
                             ##"conservation","align_score"
                             ##"profileSD.m","profileSD.mir"
                             )

if( data.organism=="simulation" ){
  ## info$interactionfactors <- c("constant",
  ##                              "pred")
  info$interactionfactors <- c("constant")
}  
                            

if( exists("prediction.weight") ){
  info$prediction.weight <- prediction.weight
}



## #########################################


## ####################################
names(info$interactionfactors) <- info$interactionfactors

if( is.null(typefollows) && !is.null(reference) ){
  temp <- list("prior")
  names(temp) <- reference
  nonreftypes <- setdiff(unique(xpr.orig$mir$targets$type),reference)
  temp2 <- as.list( rep(reference,length(nonreftypes)) )
  names(temp2) <- nonreftypes
  typefollows <- c( temp, temp2 )
}

if( any(duplicated(names(typefollows))) ){
  stop("There are duplicate names in *typefollows*.\n")
}

typefollows <- lapply(typefollows,function(typ){
  names(typ) <- typ
  return(typ)
})

types <- names(typefollows)
names(types) <- types
info$types <- types

if( all(names(typefollows) %in% c("mir","m")) ){
  info$typefollows <- typefollows
}else{
  info$typefollows$m <- typefollows
  info$typefollows$mir <- lapply( info$types,function(x)
                                 c(prior="prior") )
}
  
  
typeisbefore <- list()
typeisbefore$mir <- lapply(types,function(typ){
  types[ which( sapply(info$typefollows$mir,function(tf) typ %in% tf ) ) ]
})
typeisbefore$m <- lapply(types,function(typ){
  types[ which( sapply(info$typefollows$m,function(tf) typ %in% tf ) ) ]
})
info$typeisbefore <- typeisbefore

xpr.orig$mir <- xpr.orig$mir[,xpr.orig$mir$targets$type %in% info$types]
xpr.orig$m <- xpr.orig$m[,xpr.orig$m$targets$type %in% info$types]

genes.orig <- list(mir=unique(xpr.orig$mir$genes[[info$gname$mir]]),
                   m=unique(xpr.orig$m$genes[[info$gname$m]]))

## possmir.genes.orig <- lapply(c(mir="mir",m="m"),function(side){
##   sapply(genes.orig[[side]],function(gn){
##     ind <- which(xpr.orig[[side]]$genes[,info$gname[[side]]]==gn)
##     ret <- xpr.orig[[side]]$genes[ind[1],info$possmir.gname[[side]]]
##   })
## })


## ###########################################

## load the prediction database (targetScan, miRanda, etc)
## data (or files previously reduced to throw out
## data we don't need) and intersect gene lists

cat("Loading prediction database information.\n")
thefilename <- paste(generatedfiledir,"possmirs_orig.Rdata",sep="/")
if( data.organism=="simulation"){
  ## do nothing
}else if( file.exists(thefilename) ){
  load( thefilename )
}else{
  source(paste(sourcedir,"loadDatabasePredictions.R",sep="/"))
}


allmirs.db <- unique(unlist(lapply(possmirs.orig,function(x)
                                x[,"miRNA.main"] )))
names(allmirs.db) <- allmirs.db

allmrnas.db <- names(possmirs.orig)
names(allmrnas.db) <- allmrnas.db

##usegind <- 1:100
##usegind <- sort(unique(round(runif(50)*200)+1))


genes <- list()
if( predictedonly ){
  genes$mir <- intersect(unique(xpr.orig$mir$genes[,info$possmir.gname$mir]),
                         allmirs.db) #[usegind]
  genes$m <- intersect(unique(xpr.orig$m$genes[,info$possmir.gname$m]),
                       allmrnas.db )#[usegind]
}else{
  genes$mir <- unique(xpr.orig$mir$genes[,info$possmir.gname$mir])
  genes$m <- unique(xpr.orig$m$genes[,info$possmir.gname$m])
}

genes$mir <- setdiff(genes$mir,c(""))
genes$m <- setdiff(genes$m,c(""))

names(genes$mir) <- genes$mir
g.mir <- length(genes$mir)
names(genes$m) <- genes$m
g.m <- length(genes$m)


cat("In total, there are",g.mir,"miRs and",
    g.m,"mRNAs under consideration.\n")


## ###########################################
## get a list of probes/genes that are "present" with significance
## *pval* in at least *numtypes* types/stages

cat("Checking probes for significant difference from dark probes.\n")
if( any(xpr.orig$mir$genes[[info$gname$mir]] %in%
        c("NegativeControl","DarkCorner")) ){
  thefilename <- paste(generatedfiledir,"presentprobes.Rdata",sep="/")
  if( file.exists( thefilename ) ){
    load( thefilename )
  }else{
    presentprobes <- lapply(c(mir="mir",m="m"),function(side){
      getPresentProbes(xpr.orig,info,side,
                       pval=0.01,numtypes=3)
    })
    if( data.organism != "simulation" ){
      save(presentprobes,file=thefilename)
    }
  }
}else{
  ## just use all probes
  presentprobes <- lapply(c(mir="mir",m="m"),function(side){
    unique( xpr.orig[[side]]$genes[,info$gname[[side]]] )
  })
}


cat("There are",length(presentprobes$mir),"miR probes and",
    length(presentprobes$m),"mRNA probes present in at least 3 sample types.\n")


## ## #######################################
## ## find probes with a high signal-to-noise ratio (SNR)
## 
## ## ## usegenelist <- genes.orig #use if *presentprobes* doesn't exist
## ## usegenelist <- presentprobes
## ## signaltonoise <- getSignalToNoiseRatio(xpr.orig,info,usegenelist)
## 
## ##save(signaltonoise,file="signaltonoise.Rdata")
## load(paste(generatedfiledir,"signaltonoise.Rdata",sep="/"))
## 
## ## x11()
## ## plot(density(signaltonoise$m))
## 
## ## take only the genes above a certain SNR
## highsnrgenes <- list()
## for( side in c("mir","m")){
##   highsnrgenes[[side]] <- names(signaltonoise[[side]])[signaltonoise[[side]]>=SNRthresh]
## }
## 
## cat("There are",length(highsnrgenes$mir),"miRs and",
##     length(highsnrgenes$m),"mRNAs with signal-to-noise ratio (SNR)above",
##     SNRthresh,".\n")


## ########################################
## use an F-statistic to test for general differential expression

thefilename <- paste(generatedfiledir,"genefstat.Rdata",sep="/")
if( file.exists(thefilename) ){
  load( thefilename )
}else{  
  genefstat <- list()
  for( side in c("mir","m")){
    genefstat[[side]] <- sapply( 1:dim(xpr.orig[[side]]$E)[1],function(i){
      dat <- data.frame( expr=xpr.orig[[side]]$E[i,],
                        type=xpr.orig[[side]]$targets$type )
      ret.lm <- lm(expr ~ type,dat)
      ret.anova <- anova( ret.lm )
      ret <- ret.anova$'Pr(>F)'[1]
      if( is.na(ret) ){
          ret <- 1/var(dat$expr)
      }
      return(ret)
    })
    names(genefstat[[side]]) <- xpr.orig[[side]]$genes[,info$gname[[side]]] 
  }
  if( data.organism != "simulation" ){
    save(genefstat,file=thefilename)
  }
}

## lapply(genefstat,length)
## lapply(genefstat,function(x) sum(x<0.05))
## lapply(signaltonoise,function(x) sum(x>5))

if( length(ftestthresh)==1 ){
  tmp <- ftestthresh
  ftestthresh <- list()
  ftestthresh$mir <- tmp
  ftestthresh$m <- tmp
}

ftestthresh$mir <- min(ftestthresh$mir,g.mir)
ftestthresh$m <- min(ftestthresh$m,g.m)


ftestgenes <- list()
genefstat.FDR <- list()
for( side in c("mir","m")){
    genefstat.FDR[[side]] <- p.adjust(genefstat[[side]],method=padjustmethod)
    if( ftestthresh[[side]] < 1 ){
    ##ftestgenes[[side]] <- names(genefstat[[side]])[genefstat[[side]]<ftestthresh]
    ftestgenes[[side]] <- names(genefstat.FDR[[side]])[genefstat.FDR[[side]]<ftestthresh[[side]]]
  }else{
    ftestgenes[[side]] <- names(sort(genefstat.FDR[[side]]))[1:ftestthresh[[side]]]
  }    
}

cat("There are",length(unique(ftestgenes$mir)),"miRs (",
    length(ftestgenes$mir),"probes ) and",
    length(unique(ftestgenes$m)),"mRNAs (",
    length(ftestgenes$m),"probes )\n  with F-test p-value less than the threshold.\n")

if( length(ftestgenes$mir) > 150 || length(ftestgenes$m) > 500 ){
  cat("\nThere are a lot of miRs or genes satisfying the f-test.\n",
      " To save computational time, it might be a good idea to change the threshold.\n",
      " The 100th miR p-value is",sort(genefstat.FDR$mir)[100],
      "and the 300th mRNA p-value is",sort(genefstat.FDR$m)[300],"\n\n")
}



## ## #################################################
## load validated targets from miRecords and TarBase
mirecords <- read.csv(paste(databasefiledir,"validated/miRecords_version3.csv",sep="/"),
                      stringsAsFactors=FALSE)
tarbase <- read.csv(paste(databasefiledir,"validated/TarBase_V5.0.csv",sep="/"),
                    stringsAsFactors=FALSE)

## ## #############
## ## make a file of all miRs on the array for input into miRWalk
## mirlist <- as.matrix(unique(xpr.orig$mir$genes[,info$gname$mir]))
## write.table(mirlist,file="mirlist.txt",row.names=FALSE,
##             col.names=FALSE,quote=FALSE)

## load miRWalk validated targets
miRWalk.targets <- read.delim(paste(databasefiledir,
                                    "validated/miRWalk_validatedTargets.txt",
                                    sep="/"), row.names=NULL,stringsAsFactors=FALSE)


## #######
mrnames <- list(mir="miRNA_mature_ID",
                m="Target.gene_name")
tbnames <- list(mir="miRNA.mod",
                m="Gene")
mwnames <- list(mir="MicroRNA.Name",
                m="Gene.Name")


## modify the names to match the data's
if( data.organism=="mouse" ){
  tarbase <- tarbase[tarbase$Organism=="Mouse",]
  tarbase$miRNA.mod <- sapply(tarbase$miRNA,function(x)
                              paste("mmu-",x,sep=""))
}else if( data.organism=="human" ){
  tarbase <- tarbase[tarbase$Organism=="Human",]
  tarbase$miRNA.mod <- sapply(tarbase$miRNA,function(x)
                              paste("hsa-",x,sep=""))
}



## ########################################
## ########################################
## cut data down to only good data and/or targets

goodlist <- list()
for( side in c("mir","m")){
  goodlist[[side]] <- intersect(presentprobes[[side]],
                                ftestgenes[[side]])
}


xpr <- NULL
if( predictedonly ){
  usegind.mir <- which( ( xpr.orig$mir$genes[[info$possmir.gname$mir]] %in%
                         genes$mir) &
                       ( xpr.orig$mir$genes[[info$gname$mir]] %in%
                        goodlist$mir ) )
  xpr$mir <- xpr.orig$mir[ usegind.mir,]


  mrna.possmirs <- genes$m[ sapply(possmirs.orig[genes$m],function(xtab){
    any(xtab[,"miRNA.main"] %in% xpr$mir$genes[,info$possmir.gname$mir])
  }) ]

  usegind.m <- which( ( xpr.orig$m$genes[[info$possmir.gname$m]] %in%
                       mrna.possmirs ) &
                     ( xpr.orig$m$genes[[info$gname$m]] %in%
                      goodlist$m ) )
  xpr$m <- xpr.orig$m[ usegind.m,]

}else{
  usegind.mir <- which( ## ( xpr.orig$mir$genes[[info$possmir.gname$mir]] %in%
                       ##  genes$mir) &
                       ( xpr.orig$mir$genes[[info$gname$mir]] %in%
                        goodlist$mir ) )
  xpr$mir <- xpr.orig$mir[ usegind.mir,]
 
  
  ## mrna.possmirs <- genes$m[ sapply(possmirs.orig[genes$m],function(xtab){
  ##   any(xtab[,"miRNA.main"] %in% xpr$mir$genes[,info$possmir.gname$mir])
  ## }) ]
  
  usegind.m <- which( ## ( xpr.orig$m$genes[[info$possmir.gname$m]] %in%
                     ##  mrna.possmirs ) &
                     ( xpr.orig$m$genes[[info$gname$m]] %in%
                      goodlist$m ) )
  xpr$m <- xpr.orig$m[ usegind.m,]
}


## change the gene name list to match the info$gname field (possibly
## "ProbeName")



## ######################################################
## add gene information to parameter list *info*
info$genes <- list(mir=unique(xpr$mir$genes[[info$gname$mir]]),
                   m=unique(xpr$m$genes[[info$gname$m]]))

info$possmir.genes <- lapply(c(mir="mir",m="m"),function(side){
  sapply(info$genes[[side]],function(gn){
    ind <- which(xpr[[side]]$genes[,info$gname[[side]]]==gn)
    ret <- xpr[[side]]$genes[ind[1],info$possmir.gname[[side]]]
  })
})

cat("We are including",length(info$genes$mir),
    "total miR candidates\n and",
    length(info$genes$m),"mRNAs which might be targeted.\n")


## ########################################
## rescale expressions to standard normal
## and save mean expression
cat("Rescaling expression values...\n")
meanexpressions <- list()
genesdevs <- list()
maxexpressions <- list()
maxexpstages <- list()
for( side in c("mir","m") ){
  ##cat(side,"\n")
  meanexpressions[[side]] <- NULL
  genesdevs[[side]] <- NULL
  for( gn in info$genes[[side]] ){
    ind <- which( xpr[[side]]$genes[[info$gname[[side]]]]==gn )
    mat <- xpr[[side]]$E[ind,]
    mn <- mean(mat)

    stageexps <- sapply(info$types,function(typ){
      mean(xpr[[side]]$E[ind,xpr[[side]]$targets$type==typ])
    })
    maxexp.stage <- names(stageexps)[which.max(stageexps)]
    maxexp <- stageexps[maxexp.stage]
    
    sdev <- sd( as.vector(mat) )
    newmat <- (mat-mn)/sdev
    xpr[[side]]$E[ind,] <- newmat

    meanexpressions[[side]][gn] <- mn
    genesdevs[[side]][gn] <- sdev
    maxexpressions[[side]][gn] <- maxexp
    maxexpstages[[side]][gn] <- maxexp.stage
  }
}

## ## save(xpr,file="xpr_rescaled.Rdata")
## load("xpr_rescaled.Rdata")


## ##############################################
## compile the full *possmirs* list of targeting parameters

cat("Compiling the predicted miR-mRNA targeting parameters.\n")

possmir.genes <- lapply(c(mir="mir",m="m"),function(side){
  ret <- sapply(info$genes[[side]],function(gn){
    ind <- which(xpr.orig[[side]]$genes[,info$gname[[side]]]==gn)
    ret <- xpr.orig[[side]]$genes[ind[1],info$possmir.gname[[side]]]
  })
  return( ret[ret %in% genes[[side]]] )
})

mirFamiliesGCcontent <- read.delim(paste(databasefiledir,
                                         "predicted/targetscan/miR_Family_Info.txt",
                                         sep="/"),stringsAsFactors=FALSE)
ind <- which(mirFamiliesGCcontent$MiRBase.ID %in%
             possmir.genes$mir)
mirFamGC <- mirFamiliesGCcontent[ind,]

## usegenes <- intersect(possmir.genes$mir,
##                       mirFamGC$MiRBase.ID)
## names(usegenes) <- usegenes
## mirGCandLen <- t(sapply(usegenes,function(xmir){
##   ##cat(xmir,"\n")
##   mirseq <- mirFamGC$Mature.sequence[mirFamGC$MiRBase.ID==xmir]
##   splitseq <- strsplit(mirseq,split="")[[1]]
##   seqlen <- length(splitseq)
##   gccont <- sum( splitseq %in% c("G","C") ) / seqlen
##   miRfam <- mirFamGC$miR.family[mirFamGC$MiRBase.ID==xmir]
##   return(c( mirGCcontent=gccont,
##            miRlength=seqlen,
##            miRfam=miRfam) )
## }))

usegenes <- possmir.genes$m
names(usegenes) <- usegenes
possmirs <- lapply(usegenes,function(xgn){
  ##cat(xgn,"\n")
  xmat <- possmirs.orig[[xgn]]
  ##names(usecols) <- usecols
  ind <- which(xmat[,"miRNA.main"] %in% possmir.genes$mir)
  xmat.subset <- as.matrix(xmat[ind,,drop=FALSE])

  missingmirs <- unique(possmir.genes$mir[ !(possmir.genes$mir %in%
                                             xmat[,"miRNA.main"]) ])
  missingmirmat <- matrix(0,length(missingmirs),ncol(xmat))
  colnames(missingmirmat) <- colnames(xmat)
  missingmirmat[,"miRNA.main"] <- missingmirs
  missingmirmat[,"mRNA.main"] <- xgn

  if( nrow(xmat.subset)==0 ) xmat.subset <- NULL
  xmat.add <- rbind(xmat.subset,
                    missingmirmat)

  meanlogexp.mir <- sapply(xmat.add[,"miRNA.main"],function(xmir){
    pnames <- names(possmir.genes$mir)[possmir.genes$mir==xmir]
    return( mean(meanexpressions$mir[pnames]) )
  })
  pnames <- names(possmir.genes$m)[possmir.genes$m==xgn]
  meanlogexp.m <- rep( mean( meanexpressions$m[pnames] ), nrow(xmat.add) )

  
  profileSD.mir <- sapply(xmat.add[,"miRNA.main"],function(xmir){
    pnames <- names(possmir.genes$mir)[possmir.genes$mir==xmir]
    return( mean(genesdevs$mir[pnames]^2)^0.5 )
  })
  pnames <- names(possmir.genes$m)[possmir.genes$m==xgn]
  profileSD.m <- rep( mean(genesdevs$m[pnames]^2)^0.5, nrow(xmat.add) )

  maxlogexp.mir <- sapply(xmat.add[,"miRNA.main"],function(xmir){
    pnames <- names(possmir.genes$mir)[possmir.genes$mir==xmir]
    return( max(maxexpressions$mir[pnames]) )
  })

  maxlogexp.stage.mir <- sapply(xmat.add[,"miRNA.main"],function(xmir){
    pnames <- names(possmir.genes$mir)[possmir.genes$mir==xmir]
    return( maxexpstages$mir[pnames[1]] )
  })

  ## add in more columns
  ret <- cbind( xmat.add,
               constant=rep(1,nrow(xmat.add)),
               meanlogexp.mir=meanlogexp.mir,
               meanlogexp.m=meanlogexp.m,
               maxlogexp.mir=maxlogexp.mir,
               maxlogexp.stage.mir=maxlogexp.stage.mir,
               profileSD.mir=profileSD.mir,
               profileSD.m=profileSD.m) #,
                                        #mirGCandLen[xmat.add[,"miRNA.main"],,drop=FALSE])

  return(ret)
})


## ###########################################
## take out factor that aren't present in possmirs
info$interactionfactors <- intersect(info$interactionfactors,
                                     colnames(possmirs[[1]]) )
names(info$interactionfactors) <- info$interactionfactors


## Parting message
if( data.organism=="simulation" ){
  cat("Data prep has been completed for the simulation.\n")
}else{
  cat("There are",sum(sapply(possmirs,function(x){
    sum(as.numeric(as.matrix(x[,c("TS.pred","miRanda.pred")])))
  }) ), "predicted targeting interactions.\n")
    
}
