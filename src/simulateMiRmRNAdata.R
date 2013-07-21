require(limma)

predictedonly <- FALSE

possmir.gname$mir <- "ProbeName"
possmir.gname$m <- "ProbeName"


info <- NULL
info$nclus <- simnclus


## using stuff from user-defined parameters
typefollows <- lapply(typefollows,function(typ){
  names(typ) <- typ
  return(typ)
})

types <- names(typefollows)
names(types) <- types
info$types <- types
info$typefollows <- typefollows


simulateData <- function(info){

  genevar <- list()
  genevar$mir <- 1
  genevar$m <- 1e-5 # not including miR effect

  clusvar <- 1e-5
  slidevar <- 1e-5
  spotvar <- 0.1

  
  genesperclus <- 1
  slidespertype <- 3
  spotsperslide <- 1

  thresh.normal <- 0.9

  prob.notpredicted <- 0.9
  thresh.predicted <- 0.5

  ## mirtrend <- 0

  intermat.orig <- sapply(1:info$nclus$mir,function(nmir){
    ifelse( runif(info$nclus$m)>thresh.normal, 1,0 ) *
      -rgamma(info$nclus$m,1)
  })
  predmat <- sapply(1:info$nclus$mir,function(nmir){
    ifelse( runif(info$nclus$m)>prob.notpredicted, 1,0 )
  })
  intermat.pred <- sapply(1:info$nclus$mir,function(nmir){
    ifelse( runif(info$nclus$m)>thresh.predicted, 1,0 ) * 
    -rgamma(info$nclus$m,1)
  })

  intermat.final <- ifelse( predmat>0, intermat.pred, intermat.orig)
  ## intermat.final <- diag(-1,info$nclus$m,info$nclus$mir) # for debugging

  xpr.orig <- NULL
  for( side in c("mir","m") ){
    xpr.orig[[side]] <- NULL

    mireffect <- sapply(info$types,function(typ)
                        rep(0,info$nclus[[side]]))
    intermat <- matrix(0,info$nclus$m,info$nclus$mir)
    if( side=="m" ){
      ## intermat[1,1] <- 5
      intermat <- intermat.final
      
      for( typ in info$types ){
        if( "prior" %in% info$typefollows[[typ]] ){
          mireffect[,typ] <- intermat %*% clusprofs[,typ]
        }else{
          parenteffects <- intermat %*%
            sapply(info$typefollows[[typ]],function(ptyp){
              clusprofs[,typ] - clusprofs[,ptyp]
            })
          mireffect[,typ] <- apply(parenteffects,1,mean)
        }
      }
    }

    ## trend <- sapply(info$types,function(typ)
    ##                 rep(0,info$nclus[[side]]))
    ## for( typ in info$types ){
    ##   if( "prior" %in% info$typefollows[[typ]] ){
    ##     trend[,typ] <- rep(0,info$nclus[[side]])
    ##   }else{
    ##     if( side=="mir" ){
    ##       trend[,typ] <- mirtrend*c(1:info$nclus[[side]]) +
    ##         apply( trend[,info$typefollows[[typ]],drop=FALSE],1,mean )
    ##     }else{
    ##       trend[,typ] <- 0 + ##0.25*c(1:info$nclus[[side]]) +
    ##         apply( trend[,info$typefollows[[typ]],drop=FALSE],1,mean )
    ##     }
    ##   }
    ## }

    clusprofs <- sapply(info$types,function(typ)
                        rep(0,info$nclus[[side]]))
    for( typ in info$types ){
      if( "prior" %in% info$typefollows[[typ]] ){
        clusprofs[,typ] <- ##trend[,typ] + mireffect[,typ] +
          rnorm(info$nclus[[side]],sd=genevar[[side]]^0.5)
      }else{
        clusprofs[,typ] <- apply( clusprofs[,info$typefollows[[typ]],drop=FALSE],
                                 1,mean ) +
                                   ## trend[,typ] +
                                     mireffect[,typ] +
                                       rnorm(info$nclus[[side]],
                                             sd=genevar[[side]]^0.5)
      }
    }
  

    ## clusprofs <- 5*sapply(info$types,function(typ)
    ##                        rnorm(info$nclus[[side]],sd=genevar^0.5)) +
    ##                          trend + mireffect
    ## if( side=="mir" ){
    ##   clusprofs[,1] <- clusprofs[,1]*10
    ## }
   
    slidemeanprofs <- sapply(info$types,function(typ){
      rep(clusprofs[,typ],each=genesperclus) +
        rnorm(info$nclus[[side]]*genesperclus,sd=clusvar^0.5)
    })

    info$genes[[side]] <- paste(side,as.character(1:dim(slidemeanprofs)[1]),
                                sep="") #info$genes[[side]][1:dim(slidemeanprofs)[1]]
    rownames(slidemeanprofs) <- info$genes[[side]]

    xpr.orig[[side]]$targets <- data.frame(type=rep(info$types,each=slidespertype),
                                      stringsAsFactors=FALSE)
    spotmeanprofs <- sapply(xpr.orig[[side]]$targets$type,function(typ){
      slidemeanprofs[,typ] + rnorm(dim(slidemeanprofs)[1],sd=slidevar^0.5)
    })
    rownames(spotmeanprofs) <- rownames(slidemeanprofs)

    spotvals <- apply(spotmeanprofs,2,function(x){
      rep(x,each=spotsperslide) +
        rnorm(dim(slidemeanprofs)[1]*spotsperslide,sd=spotvar^0.5)
    })
    rownames(spotvals) <- rep(rownames(spotmeanprofs),
                              each=spotsperslide)

    xpr.orig[[side]]$E <- spotvals
    xpr.orig[[side]]$genes <- data.frame(ProbeName=rownames(spotvals),
                                    stringsAsFactors=FALSE)
  }


  ##xpr.orig$m <- xpr.orig$m[1:10,]
  xpr.orig$mir <- new("EList",xpr.orig$mir)
  xpr.orig$m <- new("EList",xpr.orig$m)

  colnames(intermat.final) <- info$genes$mir
  rownames(intermat.final) <- info$genes$m
  
  
  ## make possmirs.orig for simulated predictions
  colnames(predmat) <- info$genes$mir
  rownames(predmat) <- info$genes$m

  info$genes <- lapply(info$genes,function(x){
    names(x) <- x
    return(x)
  })

  possmirs.orig <- lapply(info$genes$m,function(gn){
    targmirs <- info$genes$mir[ predmat[gn,]!=0 ]
    if( length(targmirs)>0 ){
      retmat <- t( sapply(targmirs,function(tm){
        c( miRNA.main=as.vector(tm), mRNA.main=as.vector(gn), pred="1" )
      }) )
    }else{
      retmat <- t(c( miRNA.main=as.vector("none"),
                    mRNA.main=as.vector(gn), pred="1" ))[0,]
    }
    return(retmat)
  })
 
  
  return(list(info=info,
              xpr.orig=xpr.orig,
              intermat=intermat.final,
              possmirs.orig=possmirs.orig))
}



trapz <- function(mat){
  n <- nrow(mat)
  ret <- sum( diff(mat[,1]) * ( mat[1:(n-1),2] + mat[2:n,2] )/2 )
  return(ret)
}

checkInteractionResults <- function( inferredints, goldstd ){

  ## evaluate sensitivity and specificity
  simints <- goldstd ##getInteractionSignificance( goldstd, 0, eye(size(par.exper{1}.mumean,2)) );
  ##[ inferredints Ssigs ] <- getInteractionSignificance( par.Smean, par.Sprec, par.mship );
  
 
  infintvec <- abs( as.vector(inferredints) )
  intorder <- order( infintvec, decreasing=TRUE )

  simintvec <- as.vector(simints)
  sortsimints <- simintvec[intorder] != 0

  ##for( topi in 1:(nclus^2-nclus) ){

  ret <- t( sapply( 1:length(infintvec), function(topi){
    c( rocx= 1 - sum(!sortsimints[1:topi])/sum(!sortsimints), ## 1-FPR
      rocy= sum(sortsimints[1:topi])/sum(sortsimints), ## TPR
      prx= sum(sortsimints[1:topi])/sum(sortsimints), ## TPR
      pry= sum(sortsimints[1:topi])/topi ) ## PPV
  }) )

  ## make the curve complete at the ends
  ret <- rbind( c(1,0,0,1), ret, c(0,1,1,0) )

  ## x11()
  ## plot(x=ret[,c("rocx","prx")],
  ##      y=ret[,c("rocy","pry")],type="l")

  auroc <- -trapz(ret[,c("rocx","rocy")])
  aupr <- trapz(ret[,c("prx","pry")])
  
  return(list(auroc=auroc,
              aupr=aupr,
              inferredints=inferredints) )
}


## ########################
## actually runs the code and saves results to the workspace
simret <- simulateData(info)

xpr.orig <- simret$xpr.orig
info <- simret$info
possmirs.orig <- simret$possmirs.orig
intermat.true <- simret$intermat
