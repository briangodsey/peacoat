
getExpressedProbes <- function(xpr,chiptype,info,side=NULL){
  ## chiptype <- "GMP"
  cat("--",chiptype,"\n")
  ind <- which(## xpr[[side]]$targets$type %in% c("STHSC","LTHSC") &
               xpr[[side]]$targets$type %in% c(chiptype) ) # &
                                        # xpr[[side]]$targets$RNA=="miR" &
                                        # xpr[[side]]$targets$ITD %in% c(0))
  
  ## nd <- lapply(ind,function(x){
  ##   dat <- log(xpr[[side]]$E[,x])
  ##   x11(); hist(dat[dat>5])
  ## })
  
  dat <- xpr[[side]]$E[,ind,drop=FALSE]
  negcontrols <- xpr[[side]]$E[xpr[[side]]$genes[[info$gname[[side]]]] %in%
                               c("NegativeControl","DarkCorner"),
                       ind,drop=FALSE]
  ##negcontrols <- log(xpr[[side]]$Eb[,ind,drop=FALSE])

  ## do it namewise
  sysnames <- unique(xpr[[side]]$genes[[info$gname[[side]]]])
  names(sysnames) <- sysnames
  probp.name <- t(sapply(sysnames,function(sn){
    i <- which(xpr[[side]]$genes[[info$gname[[side]]]]==sn)
    ##cat(length(i),"-",dat[i,],sn,"\n")
    Pval <- t.test(dat[i,],negcontrols,var.equal=TRUE,
                   alternative="greater")$p.value
    meanlogexp <- mean( as.vector(dat[i,]) )
    varlogexp <- var( as.vector(dat[i,]) )
    ret <- unlist(list(Pval=Pval,meanlogexp=meanlogexp,
                       meanexp=exp(meanlogexp),varlogexp=varlogexp))
    return(ret)
  }))
  lowpvaltab.name <- probp.name[order(probp.name[,"meanexp"],
                                      decreasing=TRUE),]
  ## sum(lowpvaltab.name[,"Pval"]<0.01)
  ## write.csv(lowpvaltab.name,file="namewisePresence_GMP_withoutITD.csv",
  ##          quote=FALSE)
  return(lowpvaltab.name)
}


## get a list of probes/genes that are "present" with significance
## *pval* in at least *numtypes* types/stages
getPresentProbes <- function(xpr.orig,info,side,
                             pval=0.01,numtypes=3){
  cat(side,"\n")
  pprobes <- lapply(types,function(typ){
    tab <- getExpressedProbes(xpr.orig,typ,info,side=side)
    return( rownames(tab[tab[,"Pval"]<pval,]) )
  })

  ## ## convert from probes to genes
  ## pgenes <- lapply(pprobes,function(x){
  ##   unique(xpr.orig$m$genes$GeneName[ xpr.orig$m$genes$ProbeName %in% x ])
  ## })
  ## as.matrix(sapply(pgenes,length))
  
  numpresent <- table(unlist(pprobes))
  presentprobes <- names(numpresent)[numpresent>=numtypes]
  return(presentprobes)
}


## calculate the signal (gene profile standard deviation) to noise
## (standard deviation of slide means) ratios
getSignalToNoiseRatio <- function(xpr.orig,info,usegenelist){
  highsnrgenes <- NULL
  signaltonoise <- NULL
  for( side in c("mir","m")){
    ## allspotvar <- sapply(usegenelist[[side]],function(gn){
    ##   ind <- which(xpr.orig[[side]]$genes[[info$gname[[side]]]]==gn)
    ##   apply(xpr.orig[[side]]$E[ind,,drop=FALSE],2,var)
    ## })
    allspotmean <- sapply(usegenelist[[side]],function(gn){
      ind <- which(xpr.orig[[side]]$genes[[info$gname[[side]]]]==gn)
      apply(xpr.orig[[side]]$E[ind,,drop=FALSE],2,mean)
    })
    allslidevar <- sapply(info$types,function(typ){
      ind <- which(xpr.orig[[side]]$targets$type==typ)
      apply(allspotmean[ind,,drop=FALSE],2,var)
    })
    alltypemean <- sapply(info$types,function(typ){
      ind <- which(xpr.orig[[side]]$targets$type==typ)
      apply(allspotmean[ind,,drop=FALSE],2,mean)
    })
    allgenevar <- apply(alltypemean,1, function(x) var(x,na.rm=TRUE))
    
    meantechvar <- apply(allslidevar,1,function(x)
                         mean(x,na.rm=TRUE))
    signaltonoise[[side]] <- (allgenevar/meantechvar)^0.5
    
    ## median(allspotvar[allspotvar>0],na.rm=TRUE)
    ## median(allslidevar[allslidevar>0],na.rm=TRUE)
    ## median(allgenevar[allgenevar>0],na.rm=TRUE)
  }
  return(signaltonoise)
}


## calculate the raw correlations between miR and mRNA profiles
getAllCorrelations <- function(xpr,info,usegenelist){

  typemeans <- sapply( c(mir="mir",m="m"),function(side){
    ## allspotvar <- sapply(usegenelist[[side]],function(gn){
    ##   ind <- which(xpr[[side]]$genes[[info$gname[[side]]]]==gn)
    ##   apply(xpr[[side]]$E[ind,,drop=FALSE],2,var)
    ## })
    allspotmean <- sapply(usegenelist[[side]],function(gn){
      ind <- which(xpr[[side]]$genes[[info$gname[[side]]]]==gn)
      apply(xpr[[side]]$E[ind,,drop=FALSE],2,mean)
    })
    allslidevar <- sapply(info$types,function(typ){
      ind <- which(xpr[[side]]$targets$type==typ)
      apply(allspotmean[ind,,drop=FALSE],2,var)
    })
    alltypemean <- sapply(info$types,function(typ){
      ind <- which(xpr[[side]]$targets$type==typ)
      apply(allspotmean[ind,,drop=FALSE],2,mean)
    })
    return(alltypemean)
  })

  seqtypechanges <- lapply( c(mir="mir",m="m"),function(side){
    ret2 <- sapply(info$types,function(typ){
      if( "prior" %in% info$typefollows[typ] ){
        ret <- rep(0,dim(typemeans[[side]])[1])
      }else{         
        ind <- which(colnames(typemeans[[side]])==typ)
        parents <- which(colnames(typemeans[[side]]) %in%
                         info$typefollows$m[[typ]])
        ret <- typemeans[[side]][,ind] -
          apply(typemeans[[side]][,parents,drop=FALSE],1,mean)
      }
    })
    rownames(ret2) <- rownames(typemeans[[side]])
    return(ret2)
  })

  allcors <- list(rawcors=cor(t(typemeans$m),
                    t(typemeans$mir),use="complete.obs"),
                  seqcors=cor(t(seqtypechanges$m),
                    t(seqtypechanges$mir),use="complete.obs"))
  return(allcors)
}

