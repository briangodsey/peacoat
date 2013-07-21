generateNiceResults <- function(xpr,llh,pars,priors,info,possmirs){

  ## fit.lm <- list()
  ## for( side in c("mir","m")){
  ##   desmat <- model.matrix(~type,xpr[[side]]$targets)
  ##   fit.lm <- lmFit(xpr[[side]],design=desmat)
  ## }

  usegenelist <- info$genes
  allcors <- getAllCorrelations(xpr,info,usegenelist)


  targetList <- lapply(info$genes$mir,function(mirgene){
    ret1 <- lapply(info$genes$m,function(mgene){
      ##cat(mirgene,mgene,"\n")
      mirname <- info$possmir.genes$mir[mirgene]
      mname <- info$possmir.gene$m[mgene]
      ind <- which(possmirs[[mname]][,"miRNA.main"]==mirname)
      if( length(ind) > 0 ){
        parmat <- possmirs[[mname]][ind,,drop=FALSE]
      }else{
        parmat <- t(rep(0,dim(possmirs[[1]])[2]))
        colnames(parmat) <- colnames(possmirs[[1]])
        parmat[,"miRNA.main"] <- mirname
        parmat[,"mRNA.main"] <- mname
        ##parmat[,"miRfam"] <- mirGCandLen[mirname,"miRfam"]
      }

      parmat.single <- data.frame(miRNA.main=parmat[1,"miRNA.main"],
                                  mRNA.main=parmat[1,"mRNA.main"],
                                  TS.pred=sum(as.numeric(parmat[,"TS.pred"])),
                                  miRanda.pred=sum(as.numeric(parmat[,"miRanda.pred"])),
                                  meanlogexp.mir=parmat[1,"meanlogexp.mir"],
                                  meanlogexp.m=parmat[1,"meanlogexp.m"],
                                  maxlogexp.mir=parmat[1,"maxlogexp.mir"],
                                  maxlogexp.stage.mir=parmat[1,"maxlogexp.stage.mir"],
                                  profileSD.mir=parmat[1,"profileSD.mir"],
                                  profileSD.m=parmat[1,"profileSD.m"],
                                  stringsAsFactors=FALSE)
      
      indint <- pars$interactions$individual[mgene,mirgene]
      indint.sd <- pars$interactions$individual.prec[mgene,mirgene]^-0.5
      indint.z <- indint/indint.sd
      ## snr.mir <- signaltonoise$mir[mirgene]
      ## snr.m <- signaltonoise$m[mgene]
      ##mirtrend <- pars$mir$trend[pars$mir$clusmship[mirgene,]>0.5]
      ##mirtrend.z <- pars$mir$trend[pars$mir$clusmship[mirgene,]>0.5] /
      ##  pars$mir$trend.prec[pars$mir$clusmship[mirgene,]>0.5]^-0.5
      mirtrend <- pars$mir$trend[ which.max(pars$mir$clusmship[mirgene,]) ]
      mirtrend.z <- pars$mir$trend[ which.max(pars$mir$clusmship[mirgene,]) ] /
        pars$mir$trend.prec[ which.max(pars$mir$clusmship[mirgene,]) ]^-0.5

      valid.mirecords <- ifelse( any(mirecords[ mirecords[,mrnames$mir]==mirname,
                                               mrnames$m ]==mname), 1, 0)

      ## if( data.organism=="mouse" ){
      ##   mirgene.add <- paste("mmu-",mirgene,sep="")
      ## }else if( data.organism=="human" ){
      ##   mirgene.add <- paste("hsa-",mirgene,sep="")
      ## }

      if( all(unlist(tbnames) %in% colnames(tarbase)) ){
        valid.tarbase <- ifelse( any(tarbase[ tarbase[,tbnames$mir]==mirname,
                                             tbnames$m ]==mname), 1, 0)
      }else{
        valid.tarbase <- 0
      }

      if( all(unlist(mwnames) %in% colnames(miRWalk.targets)) ){
        valid.miRWalk <- ifelse( any(miRWalk.targets[ miRWalk.targets[,mwnames$mir]==mirname,
                                                     mwnames$m ]==mname), 1, 0)
      }else{
        valid.miRWalk <- 0
      }

      ind <- which( miRWalk.targets[,mwnames$mir]==mirname &
                   miRWalk.targets[,mwnames$m]==mname )
      if( length(ind)>0 ){
        miRWalk.cite <- as.character(miRWalk.targets[ind[1],"Pubmed.ID"])
      }else{
        miRWalk.cite <- "NA"
      }


      ##miRfam <- ifelse(mirname %in% rownames(mirGCandLen),
      ##                 mirGCandLen[mirname,"miRfam"],"none specified")
      
      ##intmat <- data.frame(probe.mir=rep(mirgene,dim(parmat)[1]),
      ##                     probe.m=rep(mgene,dim(parmat)[1]),
      ##                     interaction=rep(indint,dim(parmat)[1]),
      ##                     interaction.SD=rep(indint.sd,dim(parmat)[1]),
      ##                     interaction.z=rep(indint.z,dim(parmat)[1]),
      ##                     miRtrend=rep(mirtrend,dim(parmat)[1]),
      ##                     miRtrend.z=rep(mirtrend.z,dim(parmat)[1]),
      ##                     rawCor=rep(allcors$rawcors[mgene,mirgene],dim(parmat)[1]),
      ##                     seqCor=rep(allcors$seqcors[mgene,mirgene],dim(parmat)[1]),
      ##                     ## snr.mir=rep(snr.mir,dim(parmat)[1]),
      ##                     ## snr.m=rep(snr.m,dim(parmat)[1]),
      ##                     valid.mirecords=rep(valid.mirecords,dim(parmat)[1]),
      ##                     valid.tarbase=rep(valid.tarbase,dim(parmat)[1]),
      ##                     valid.miRWalk=rep(valid.miRWalk,dim(parmat)[1]),
      ##                     PubMed.cite.MW=rep(miRWalk.cite,dim(parmat)[1]),
      ##                     stringsAsFactors=FALSE)
      ##fullmat <- cbind(parmat,intmat,
      ##                 stringsAsFactors=FALSE)                           

      intmat <- data.frame(probe.mir=mirgene,
                           probe.m=mgene,
                           ##miRfam=miRfam,
                           interaction=indint,
                           interaction.SD=indint.sd,
                           interaction.z=indint.z,
                           miRtrend=mirtrend,
                           miRtrend.z=mirtrend.z,
                           rawCor=allcors$rawcors[mgene,mirgene],
                           seqCor=allcors$seqcors[mgene,mirgene],
                           ## snr.mir=snr.mir,
                           ## snr.m=snr.m,
                           valid.mirecords=valid.mirecords,
                           valid.tarbase=valid.tarbase,
                           valid.miRWalk=valid.miRWalk,
                           PubMed.cite.MW=miRWalk.cite,
                           stringsAsFactors=FALSE)
      fullmat <- cbind(parmat.single,intmat,
                       stringsAsFactors=FALSE)                           

      return(fullmat)
    })
    ret2 <- ret1[[1]]
    for( i in 2:length(ret1) ){
      ret2 <- rbind(ret2,ret1[[i]])
    }
    rownames(ret2) <- NULL
    return(ret2)
  })

  targetTable <- targetList[[1]]
  for( i in 2:length(targetList) ){
    targetTable <- rbind(targetTable,targetList[[i]])
  }
  rownames(targetTable) <- NULL

  targetTable.sort.intz <- targetTable[rev(order(abs(targetTable$interaction.z))),]
  rownames(targetTable.sort.intz) <- NULL
  targetTable.sort.cor <- targetTable[rev(order(abs(targetTable$rawCor))),]
  rownames(targetTable.sort.cor) <- NULL
  targetTable.sort.seqcor <- targetTable[rev(order(abs(targetTable$seqCor))),]
  rownames(targetTable.sort.seqcor) <- NULL

  targetTable.sort.intz.neg <- targetTable[order(targetTable$interaction.z),]
  rownames(targetTable.sort.intz.neg) <- NULL
  targetTable.sort.cor.neg <- targetTable[order(targetTable$rawCor),]
  rownames(targetTable.sort.cor.neg) <- NULL
  targetTable.sort.seqcor.neg <- targetTable[order(targetTable$seqCor),]
  rownames(targetTable.sort.seqcor.neg) <- NULL

  targetTable.sort.trendz <- targetTable[rev(order(abs(targetTable$miRtrend.z))),]
  rownames(targetTable.sort.trendz) <- NULL

  return( list(raw=targetTable,
               intz=targetTable.sort.intz,
               intz.neg=targetTable.sort.intz.neg,
               cor=targetTable.sort.cor,
               cor.neg=targetTable.sort.cor.neg,
               seqcor=targetTable.sort.seqcor,
               seqcor.neg=targetTable.sort.seqcor.neg,
               trendz=targetTable.sort.trendz) )
}
