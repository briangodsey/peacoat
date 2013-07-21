## library(mcmc) ## only for alternative version of
## optimizeDevelopments()

library(MASS) # for ginv()

## #############################################
## set up initial parameters/priors

initpriors <- function(info,possmirs){
  priors <- list(
                 interactions=
                 list(## cluster=list(mean=rep(0,info$nclus$mir),
                      ##   prec=diag(rep(1e-1,info$nclus$mir))),
                      intercoeffs=list(mean=sapply(info$interactionfactors,function(x) 0),
                        prec=1e-5*diag(rep(1,length(info$interactionfactors)))),
                      ## individualinteractions=sapply(info$genes$mir,function(x)
                      ##   sapply(info$genes$m,function(y) 0 ) ),
                      ## individualinteractions.prec=sapply(info$genes$mir,function(x)
                      ##   sapply(info$genes$m,function(y) 1e5 ) ),
                      ## prec.individual=list(alf=1e-5,bta=1e-5),
                      ## prec.cluster=list(alf=1e-5,bta=1e-5),
                      prec.common=list(alf=1e-5,bta=1e-5)
                      ##prediction.var=0.1
                      ),

                 mir=
                 list(techprec=list(spot=list(alf=1e-5,bta=1e-5),
                        slide=list(alf=1e-5,bta=1e-5)),
                      clusprec=list(alf=1), #bta=1e-0),
                      clusprecBTAhyp=list(alf=10,bta=1e10),
                      stageprec=list(alf=1), #bta=1e-1),
                      stageprecBTAhyp=list(alf=1,bta=1),
                      ##devprec=list(alf=1e-5,bta=1e-5),
                      development=list(alf=1,bta=1),
                      stagemean=list(mean=0,prec=1e1),
                      clusmship=sapply(1:info$nclus$mir,function(i){
                        sapply(info$genes$mir,function(gn) 1/info$nclus$mir)
                      }),
                      trend=list(mean=0,prec=1e10)
                      ),

                 m=
                 list(techprec=list(spot=list(alf=1e-5,bta=1e-5),
                        slide=list(alf=1e-5,bta=1e-5)),
                      clusprec=list(alf=1), #bta=1e-5),
                      clusprecBTAhyp=list(alf=10,bta=1e10),
                      stageprec=list(alf=1), #bta=1e-1),
                      stageprecBTAhyp=list(alf=1,bta=1),
                      ##devprec=list(alf=1e-5,bta=1e-5),
                      development=list(alf=1,bta=1),
                      stagemean=list(mean=0,prec=1e1),
                      clusmship=sapply(1:info$nclus$m,function(typ){
                        sapply(info$genes$m,function(gn) 1/info$nclus$m)
                      }),
                      trend=list(mean=0,prec=1e10)
                      )
                 )
  return(priors)
}            


random.indicator.mat <- function(k,nams){
  n <- length(nams)
  reorderedn <- order(runif(n))
  reorderedk <- c( order(runif(k)), ceiling( runif(n-k)*k) )
  mat <- matrix(0,n,k)
  for( i in 1:n ){
    mat[reorderedn[i],reorderedk[reorderedn[i]]] <- 1
  }
  rownames(mat) <- nams
  return(mat)
}

initpars <- function(info,possmirs){
  nonpriortypes <- list()
  nonpriortypes$mir <- sapply(info$typefollows$mir,
                              function(x) !("prior" %in% x))
  nonpriortypes$m <- sapply(info$typefollows$m,
                            function(x) !("prior" %in% x))
  pars <- list(
               interactions=
               list(cluster=matrix(0,info$nclus$m,info$nclus$mir),
                    cluster.prec=lapply(1:info$nclus$m,function(xi)
                      diag( rep(1,info$nclus$mir) ) ),
                    intercoeffs=sapply(info$interactionfactors,function(x) 0),
                    intercoeffs.prec=1e10*diag(rep(1,length(info$interactionfactors))),
                    individual=sapply(info$genes$mir,function(x)
                      sapply(info$genes$m,function(y) 0 ) ),
                    individual.prec=sapply(info$genes$mir,function(x)
                      sapply(info$genes$m,function(y) 1e5 ) ),
                    ##prec.individual=list(alf=1e-5,bta=1),
                    ## prec.individual=list(alf=sapply(info$genes$mir,function(x)
                    ##                        sapply(info$genes$m,function(y) 1e-5 )),
                    ##   bta=sapply(info$genes$mir,function(x)
                    ##     sapply(info$genes$m,function(y) 1 ))),
                    ## prec.cluster=list(alf=rep(1e-5,info$nclus$mir),
                    ##   bta=rep(1,info$nclus$mir))
                    prec.common=list(alf=1e-5,bta=1)
                    ),

               mir=
               list(techprec=list(spot=list(alf=1,bta=1),
                      slide=list(alf=0.1,bta=1)),
                    clusprec=list(alf=sapply(info$types,function(typ){
                      rep(1e-3,info$nclus$mir)
                    }),
                      bta=sapply(info$types,function(typ){
                        rep(1,info$nclus$mir)
                      })),
                    clusprecBTAhyp=list(alf=1e1,bta=1e2),
                    stageprec=list(alf=rep(1e-2,info$nclus$mir),
                      bta=rep(1,info$nclus$mir)),
                    stageprecBTAhyp=list(alf=1e1,bta=1e2),
                    ##devprec=list(alf=0.1,bta=1),
                    development=lapply(info$typefollows$mir[nonpriortypes$mir],
                      function(typ) sapply(typ,function(subtyp) 1 )),
                    ##list(alf=1, #rgamma(1,10,rate=10),
                    ##   bta=1))), #rgamma(1,10,rate=10)))),
                    stagemean=sapply(info$types,function(typ){
                      rnorm(info$nclus$mir,0,sd(as.vector(xpr$mir$E)))
                    }),
                    stagemean.prec=sapply(info$types,function(typ){
                      rep(1e-1,info$nclus$mir)
                    }),
                    clusmship=random.indicator.mat(info$nclus$mir,
                      info$genes$mir),
                    ## sapply(1:info$nclus$mir,function(i){
                    ##   sapply(info$genes$mir,function(gn) runif(1)/info$nclus$mir)
                    ## }),
                    slidemean=sapply(info$types,function(typ){
                      sapply(info$genes$mir,function(gn) rnorm(1) )
                    }),
                    slidemean.prec=sapply(info$types,function(typ){
                      sapply(info$genes$mir,function(gn) 1e-1 )
                    }),
                    spotmean=apply(xpr$mir$targets,1,function(chipx){
                      sapply(info$genes$mir,function(gn) rnorm(1) )
                    }),
                    spotmean.prec=apply(xpr$mir$targets,1,function(chipx){
                      sapply(info$genes$mir,function(gn) 1e-1 )
                    }),
                    trend=rep(0,info$nclus$mir),
                    trend.prec=rep(1e10,info$nclus$mir)
                    ),

               m=
               list(techprec=list(spot=list(alf=1,bta=1),
                      slide=list(alf=0.1,bta=1)),
                    clusprec=list(alf=sapply(info$types,function(typ){
                      rep(1e-3,info$nclus$m)
                    }),
                      bta=sapply(info$types,function(typ){
                        rep(1,info$nclus$m)
                      }) ),
                    clusprecBTAhyp=list(alf=1e1,bta=1e2),
                    stageprec=list(alf=rep(1e-2,info$nclus$m),
                      bta=rep(1,info$nclus$m)),
                    stageprecBTAhyp=list(alf=1e1,bta=1e2),
                    ## devprec=list(alf=0.1,bta=1),
                    development=lapply(info$typefollows$m[nonpriortypes$m],
                      function(typ) sapply(typ,function(subtyp) 1 )),
                    ##list(alf=1, #rgamma(1,10,rate=10),
                    ##bta=1))), #rgamma(1,10,rate=10)))),
                    stagemean=sapply(info$types,function(typ){
                      rnorm(info$nclus$m,0,sd(as.vector(xpr$m$E)))
                    }),
                    stagemean.prec=sapply(info$types,function(typ){
                      rep(1e-1,info$nclus$m)
                    }),
                    clusmship=random.indicator.mat(info$nclus$m,
                      info$genes$m),
                    ## clusmship=sapply(1:info$nclus$m,function(i){
                    ##   sapply(info$genes$m,function(gn) runif(1)/info$nclus$m)
                    ## }),
                    slidemean=sapply(info$types,function(typ){
                      sapply(info$genes$m,function(gn) rnorm(1) )
                    }),
                    slidemean.prec=sapply(info$types,function(typ){
                      sapply(info$genes$m,function(gn) 1e-1 )
                    }),
                    spotmean=apply(xpr$m$targets,1,function(chipx){
                      sapply(info$genes$m,function(gn) rnorm(1) )
                    }),
                    spotmean.prec=apply(xpr$m$targets,1,function(chipx){
                      sapply(info$genes$m,function(gn) 1e-1 )
                    }),
                    trend=rep(0,info$nclus$m), #rnorm(info$nclus$mir)/10,
                    trend.prec=rep(1e10,info$nclus$m) #rgamma(info$nclus$mir,1)
                    )
               )
  return(pars)
}  


## #############################
## calculates the expectaions of gamma-distributed variables so that
## we don't waste time calculating them many times over

gamExpLogx <- function( shape, scale ){
  ent <- psigamma(shape) - log(scale)
  return(ent)
}

getExpectations <- function(pars){
  expect <- list(
                 interactions=
                 list(
                      ## prec.individual=pars$interactions$prec.individual$alf/pars$interactions$prec.individual$bta,
                      ## prec.cluster=pars$interactions$prec.cluster$alf/pars$interactions$prec.cluster$bta
                      prec.common=pars$interactions$prec.common$alf/pars$interactions$prec.common$bta
                      ),
                 mir=
                 list(techprec=list(spot=pars$mir$techprec$spot$alf/pars$mir$techprec$spot$bta,
                        slide=pars$mir$techprec$slide$alf/pars$mir$techprec$slide$bta),
                      clusprec=pars$mir$clusprec$alf/pars$mir$clusprec$bta,
                      clusprecBTAhyp=pars$mir$clusprecBTA$alf/pars$mir$clusprecBTA$bta,
                      stageprec=pars$mir$stageprec$alf/pars$mir$stageprec$bta,
                      stageprecBTAhyp=pars$mir$stageprecBTA$alf/pars$mir$stageprecBTA$bta #,
                      ## development=lapply(pars$mir$development,function(x){
                      ##   sapply(x,function(y) y$alf/y$bta )
                      ## }),
                      ## devprec=pars$mir$devprec$alf/pars$mir$devprec$bta
                      ),
                 m=
                 list(techprec=list(spot=pars$m$techprec$spot$alf/pars$m$techprec$spot$bta,
                        slide=pars$m$techprec$slide$alf/pars$m$techprec$slide$bta),
                      clusprec=pars$m$clusprec$alf/pars$m$clusprec$bta,
                      clusprecBTAhyp=pars$m$clusprecBTA$alf/pars$m$clusprecBTA$bta,
                      stageprec=pars$m$stageprec$alf/pars$m$stageprec$bta,
                      stageprecBTAhyp=pars$m$stageprecBTA$alf/pars$m$stageprecBTA$bta #,
                      ## development=lapply(pars$m$development,function(x){
                      ##   sapply(x,function(y) y$alf/y$bta )
                      ## }),
                      ## devprec=pars$m$devprec$alf/pars$m$devprec$bta
                      )
                 )
  return(expect)
}


## ################################################
## update functions

updateWithinSlide <- function( xpr, pars, priors, info, side=NULL ){
  expect <- getExpectations( pars )
  
  pars[[side]]$spotmean.prec <-t(sapply(info$genes[[side]],function(gn){
    gnxpr <- xpr[[side]]$E[xpr[[side]]$genes[[info$gname[[side]]]]==gn,,drop=FALSE]
    totwts <- expect[[side]]$techprec$slide *
      rep(1,dim(xpr[[side]]$E)[2]) +
        apply( expect[[side]]$techprec$spot *
              matrix(1,dim(gnxpr)[1],dim(gnxpr)[2]),2,sum)
    return(totwts)
  }))
  colnames(pars[[side]]$spotmean.prec) <- colnames(xpr[[side]]$E)

  wtedsums <-t(sapply(info$genes[[side]],function(gn){
    gnxpr <- xpr[[side]]$E[xpr[[side]]$genes[[info$gname[[side]]]]==gn,,drop=FALSE]
    wsum <- expect[[side]]$techprec$slide *
      pars[[side]]$slidemean[gn,xpr[[side]]$targets$type] +
        apply(expect[[side]]$techprec$spot*gnxpr,2,sum)
    return(wsum)
  }))
  pars[[side]]$spotmean <- wtedsums / pars[[side]]$spotmean.prec
  colnames(pars[[side]]$spotmean) <- colnames(xpr[[side]]$E)

  
  expect <- getExpectations( pars )

  ## ##############
  ## techprec$spot
  pars[[side]]$techprec$spot$alf <- priors[[side]]$techprec$spot$alf +
    0.5*prod(dim(xpr[[side]]$E))
  pars[[side]]$techprec$spot$bta <- priors[[side]]$techprec$spot$bta +
    0.5*sum( sapply(info$genes[[side]],function(gn){
      ind <- which(xpr[[side]]$genes[[info$gname[[side]]]]==gn)
      gsum <- sum( sapply(ind,function(xi){
        rsum <- sum( xpr[[side]]$E[xi,]^2 +
                    pars[[side]]$spotmean[gn,]^2 +
                    1/pars[[side]]$spotmean.prec[gn,] -
                    2*xpr[[side]]$E[xi,]*pars[[side]]$spotmean[gn,] )
      }) )
    }) )
  
  return(pars)
}


updateAcrossSlides <- function( xpr, pars, priors, info, side=NULL ){
  expect <- getExpectations( pars )

  pars[[side]]$slidemean.prec <- pars[[side]]$clusmship %*% expect[[side]]$clusprec +
    sapply(info$types,function(typ){
      ds <- dim(pars[[side]]$spotmean[,xpr[[side]]$targets$type==typ,drop=FALSE])
      wts <- matrix(ds[2]*expect[[side]]$techprec$slide,nrow=ds[1],ncol=1)
      return(wts)
    })
  
  wtedsums <- pars[[side]]$clusmship %*% (expect[[side]]$clusprec * pars[[side]]$stagemean) +
    sapply(info$types,function(typ){
      mat <- pars[[side]]$spotmean[,xpr[[side]]$targets$type==typ,drop=FALSE]
      wsum <- apply(mat*expect[[side]]$techprec$slide,1,sum)
      return(wsum)
    })
  pars[[side]]$slidemean <- wtedsums/pars[[side]]$slidemean.prec


  expect <- getExpectations( pars )
  
  ## ################
  ## techprec$slide
  pars[[side]]$techprec$slide$alf <- priors[[side]]$techprec$slide$alf +
    0.5 * prod( dim(pars[[side]]$spotmean) )
  pars[[side]]$techprec$slide$bta <- priors[[side]]$techprec$slide$bta +
    0.5 * sum( pars[[side]]$spotmean^2 + 1/pars[[side]]$spotmean.prec +
              pars[[side]]$slidemean[,xpr[[side]]$targets$type]^2 +
              1/pars[[side]]$slidemean.prec[,xpr[[side]]$targets$type] -
              2 * pars[[side]]$spotmean*pars[[side]]$slidemean[,xpr[[side]]$targets$type] )

  return(pars)
}


getMiReffect <- function(pars,info,side=NULL){

  if( side=="m" ){
    mireffect.mat <- pars$interactions$cluster %*% pars$mir$stagemean
  }else{
    mireffect.mat <- sapply(info$types,function(typ)
                            rep(0,info$nclus$mir) )
  }
  
  mireffect.list <- lapply(info$types,function(typ){
    prevtyp <- info$typefollows$m[[typ]]
    if( side=="m" ){
      if( "prior" %in% prevtyp[[1]] ){
        mireff <- matrix(0,info$nclus[[side]],1) #mireffect.mat[,typ,drop=FALSE]
      }else{
        mireff <- sapply( prevtyp,function(xptyp)
                         mireffect.mat[,typ] - mireffect.mat[,xptyp] )
      }
    }else{
      mireff <- sapply( prevtyp,function(xptyp)
                       rep(0,info$nclus[[side]]) ) # matrix(0,info$nclus[[side]],1)
    }
    return(mireff)
  })

  return(mireffect.list)
}

getMiReffectVAR <- function(pars,info,side=NULL,
                            intClusPrecInv=NULL ){
  expect <- getExpectations( pars )

  ## a bit of pre-calculation
  if( side=="m" && is.null(intClusPrecInv) ){
    intClusPrecInv <- lapply(pars$interactions$cluster.prec,ginv)
  }

  mireffectVAR.list <- lapply(info$types,function(typ){
    prevtyp <- info$typefollows$m[[typ]]
    if( side=="m" ){
      if( "prior" %in% prevtyp[[1]] ){
        mireffVAR <- rep(0,info$nclus[[side]]) 
      }else{
        mirdiff <- sapply( prevtyp,function(xptyp)
                          pars$mir$stagemean[,typ] -
                          pars$mir$stagemean[,xptyp] )
        mirdiffVAR <- sapply( prevtyp,function(xptyp)
                             pars$mir$stagemean.prec[,typ]^-1 +
                             pars$mir$stagemean.prec[,xptyp]^-1 )
        avgmirdiff <- as.matrix( apply(mirdiff,1,function(x){
          sum( x * abs(pars$m$development[[typ]])^-1 ) /
            sum(abs(pars$m$development[[typ]])^-1)
        }) )
        avgmirdiffVAR <- diag( apply(mirdiffVAR,1,function(x){
          sum( x * abs(pars$m$development[[typ]])^-2 ) /
            sum(abs(pars$m$development[[typ]])^-1)^2
        }) )
        middlemat <- avgmirdiff %*% t(avgmirdiff) + avgmirdiffVAR

        ## {Rprof()
        ##  mireffVAR1 <- sapply(1:info$nclus$m,function(xi){
        ##    ret <- pars$interactions$cluster[xi,,drop=FALSE] %*%
        ##      avgmirdiffVAR %*% ## middlemat %*%
        ##        t( pars$interactions$cluster[xi,,drop=FALSE] ) +
        ##          sum( diag( intClusPrecInv[[xi]] %*%
        ##                    middlemat ) )
        ## return(ret)
        ##  })
        ##  Rprof(NULL)}

        mireffVAR <- sapply(1:info$nclus$m,function(xi){
          ret <- pars$interactions$cluster[xi,,drop=FALSE] %*%
            avgmirdiffVAR %*% ## middlemat %*%
              t( pars$interactions$cluster[xi,,drop=FALSE] ) +
                sum( intClusPrecInv[[xi]] * t(middlemat) ) ## this is the trace of the product
          return(ret)
        })
      }
    }else{
      mireffVAR <- rep(0,info$nclus[[side]])
    }
    return(mireffVAR)
  })
  
  return(mireffectVAR.list)
}



updateStages <- function( xpr, pars, priors, info, side=NULL,
                         intClusPrecInv=NULL){

  ## llhmat <- matrix(0,info$nclus[[side]],length(info$types))
  ## colnames(llhmat) <- info$types
  ## llhmat2 <- llhmat
  ## llhmat3 <- llhmat
  ## 
  ## llhold <- calcllh(priors,pars,info,possmirs)

  ## pre-calculating
  mireffect.list <- getMiReffect(pars,info,side=side)

  if( side=="mir" && is.null(intClusPrecInv) ){
    intClusPrecInv <- lapply(pars$interactions$cluster.prec,ginv)
  }

  if( side=="mir" && !is.null(intClusPrecInv) ){
    ## this makes things faster later
    icpi.array <- array(dim=c(info$nclus$m,info$nclus$mir,info$nclus$mir))
    for( yi in 1:length(intClusPrecInv) ){
      icpi.array[yi,,] <- intClusPrecInv[[yi]]
    }
  }
  
  for( dotypes in info$types )
    {
      for( iclus in 1:info$nclus[[side]] )
        {
          expect <- getExpectations( pars )

          ## ##############
          ## precisions
          priorpart.prec <- sapply(dotypes,function(typ){
            prevtyp <- info$typefollows[[side]][[typ]]
            if( prevtyp[[1]]=="prior" ){
              precs <- priors[[side]]$stagemean$prec
            }else{
              ## the weighted average of weights from preceding stages
              devels <- pars[[side]]$development[[typ]]
              avgWT <- sum( abs(devels)^-2 ) /
                sum( abs(devels)^-1 )
              precs <- expect[[side]]$stageprec[iclus] * avgWT
            }
            return(precs)
          })

          postpart.prec <- sapply(dotypes,function(typ){
            posttyp <- info$typeisbefore[[side]][[typ]]
            if( length(posttyp)==0 ){
              precs <- 0
            }else{
              postdevWTS <- sapply( posttyp,function(ptyp){
                devels <- pars[[side]]$development[[ptyp]]
                avgWT <- sum( abs(devels)^-2 ) /
                  sum( abs(devels)^-1 )
                parenthood <- abs(devels[typ])^-1 /
                  sum( abs(devels)^-1 )
                return( avgWT * parenthood^2 )
              })
              precs <- expect[[side]]$stageprec[iclus] * sum(postdevWTS)
            }
            return(precs)
          })

          cluspart.prec <- sapply(dotypes,function(typ){
            sum( pars[[side]]$clusmship[,iclus]*
                expect[[side]]$clusprec[iclus,typ] )
          })

          interpart.prec <- sapply(dotypes,function(typ){
            if( side=="mir" ){
              prevtyp <- info$typefollows$m[[typ]]
              if( prevtyp[[1]]=="prior" ){
                precsum.aschild <- 0
              }else{
                devels <- pars$m$development[[typ]]
                wtdAvgInvDev <- sum( abs(devels)^-2 ) /
                  sum( abs(devels)^-1 )
                ## allprecs.aschild <- sapply( 1:info$nclus$m,function(yi){
                ##   wtdAvgInvDev * expect$m$stageprec[yi] *
                ##     ( pars$interactions$cluster[yi,iclus]^2  +
                ##      intClusPrecInv[[yi]][iclus,iclus] )
                ## })
                allprecs.aschild <- wtdAvgInvDev * expect$m$stageprec *
                  ( pars$interactions$cluster[,iclus]^2  +
                   icpi.array[,iclus,iclus] )
                precsum.aschild <- sum(allprecs.aschild)
              }

              posttyp <- info$typeisbefore$m[[typ]]
              if( length(posttyp)==0 ){
                precsum.asparent <- 0
              }else{
                avgPostWts <- sapply( posttyp,function(ptyp){
                  sum( abs(pars$m$development[[ptyp]])^-2) /
                    sum( abs(pars$m$development[[ptyp]])^-1 )
                })

                ## allprecs.asparent <- sapply( 1:info$nclus$m,function(yi){
                ##   ret2 <- sapply( posttyp,function(ptyp){
                ##     ret <- avgPostWts[ptyp] * expect$m$stageprec[yi] *
                ##       ( pars$interactions$cluster[yi,iclus]^2 +
                ##        intClusPrecInv[[yi]][iclus,iclus] ) *
                ##          abs(pars$m$development[[ptyp]][typ])^-2 /
                ##            sum( abs(pars$m$development[[ptyp]])^-1 )^2
                ##     return(ret)
                ##   })
                ##   return(ret2)
                ## })

                allprecs.asparent <- sapply( posttyp,function(ptyp){
                  ret <- avgPostWts[ptyp] * expect$m$stageprec *
                    ( pars$interactions$cluster[,iclus]^2 +
                     icpi.array[,iclus,iclus] ) *
                       abs(pars$m$development[[ptyp]][typ])^-2 /
                         sum( abs(pars$m$development[[ptyp]])^-1 )^2
                  return(ret)
                })
                
                precsum.asparent <- sum(allprecs.asparent)
              }
              
              precs <- precsum.aschild + precsum.asparent
            }else{
              precs <- 0
            }
            return(precs)
          })

          
          ## #########################
          ## adding in the means and multiplying by the above precisions
          
          priorpart <- sapply(dotypes,function(typ){
            prevtyp <- info$typefollows[[side]][[typ]]
            if( "prior" %in% prevtyp[[1]] ){
              prevmeans <- priors[[side]]$stagemean$mean
              precs <- priors[[side]]$stagemean$prec
            }else{
              ## the weighted average of weights from preceding stages
              devels <- pars[[side]]$development[[typ]]
              avgWT <- sum( abs(devels)^-2 ) /
                sum( abs(devels)^-1 )
              
              ## the weighted average of preceding stages by inverse
              ## development
              prevmeans <- sum( ( pars[[side]]$stagemean[iclus,unlist(prevtyp)] +
                                 mireffect.list[[typ]][iclus,] +
                                 pars[[side]]$trend[iclus] * devels ) *
                               abs(devels)^-1 ) /
                                 sum(abs(devels)^-1)
              precs <- expect[[side]]$stageprec[iclus] * avgWT
            }
            ret <- prevmeans * precs
            return(ret)
          })
          
          postpart <- sapply(dotypes,function(typ){
            posttyp <- info$typeisbefore[[side]][[typ]]
            if( length(posttyp)==0 ){
              ret <- 0
            }else{
              totsum <- sapply( posttyp,function(ptyp){
                devels <- pars[[side]]$development[[ptyp]]
                avgWT <- sum( abs(devels)^-2 ) /
                  sum( abs(devels)^-1 )
                parenthood <- abs(devels[typ])^-1 /
                  sum( abs(devels)^-1 )
                precs <- avgWT * expect[[side]]$stageprec[iclus] *
                  parenthood^2
                mireffect.self <- mireffect.list[[ptyp]][iclus,typ]
                stemeffect.self <- pars[[side]]$trend[iclus] * devels[typ]
                selfsum <- mireffect.self + stemeffect.self
                coparents <- info$typefollows[[side]][[ptyp]][which(info$typefollows[[side]][[ptyp]]!=typ)]
                if( length(coparents)>0 ){
                  mireffect.coparents <- sapply( coparents,function(cpt)
                                                mireffect.list[[ptyp]][iclus,cpt] )
                  stemeffect.coparents <- sapply( coparents,function(cpt)
                                                 pars[[side]]$trend[iclus] *
                                                 devels[cpt] )
                  coparent.wtdsum <- sum( ( pars[[side]]$stagemean[iclus,coparents] +
                                           mireffect.coparents +
                                           stemeffect.coparents ) *
                                         abs(devels)[coparents]^-1 /
                                         abs(devels)[typ]^-1 )
                  
                }else{
                  coparent.wtdsum <- 0
                }
                posttypsum <- precs *
                  ( pars[[side]]$stagemean[iclus,ptyp]/parenthood -
                   (selfsum + coparent.wtdsum) )
                return(posttypsum)
              })
              ret <- sum( totsum )
            }
            return( ret )
          })
          
          cluspart <- sapply(dotypes,function(typ){
            sum(pars[[side]]$clusmship[,iclus] * pars[[side]]$slidemean[,typ] *
                expect[[side]]$clusprec[iclus,typ])
          })

          interpart <- sapply(dotypes,function(typ){
            if( side=="mir" ){
              prevtyp <- info$typefollows$m[[typ]]
              if( prevtyp[[1]]=="prior" ){
                meansum.aschild <- 0
              }else{
                devels <- pars$m$development[[typ]]
                wtdAvgInvDev <- sum(abs(devels)^-2) /
                  sum( abs(devels)^-1 )
                ## otherstagemeans <- pars$mir$stagemean[,typ]
                ## otherstagemeans[iclus] <- 0 # we're estimating this value

                ## ## old and SLOW
                ## {Rprof()
                ## allmeans.aschild <- sapply( 1:info$nclus$m,function(yi){
                ##   wtdAvgInvDev * expect$m$stageprec[yi] *
                ##     ( pars$interactions$cluster[yi,iclus] *
                ##      ( pars$m$stagemean[yi,typ] -                      
                ##       sum( sapply(prevtyp,function(ptyp){
                ##         abs(devels[ptyp])^-1 *
                ##           ( pars$m$stagemean[yi,ptyp] +
                ##            pars$m$trend[yi]*devels[ptyp] ) /
                ##              sum(abs(devels)^-1)
                ##       }) ) ) -
                ##      sum( sapply( c(1:info$nclus$mir)[-iclus],function(omir){
                ##        ( pars$interactions$cluster[yi,iclus] *
                ##         pars$interactions$cluster[yi,omir] +
                ##         intClusPrecInv[[yi]][iclus,omir] ) *
                ##           pars$mir$stagemean[omir,typ]
                ##      }) ) - 
                ##      sum( sapply(prevtyp,function(ptyp){
                ##        sum( ( pars$interactions$cluster[yi,iclus] *
                ##              pars$interactions$cluster[yi,] +
                ##              intClusPrecInv[[yi]][iclus,] ) *
                ##            (-pars$mir$stagemean[,ptyp]) ) *
                ##              abs(devels[ptyp])^-1 / sum(abs(devels)^-1)
                ##      }) ) )
                ## })
                ##  Rprof(NULL)}

                omirind <- c(1:info$nclus$mir)[-iclus]

                ## allmeans.aschild1 <- sapply( 1:info$nclus$m,function(yi){
                ##   wtdAvgInvDev * expect$m$stageprec[yi] *
                ##     ( pars$interactions$cluster[yi,iclus] *
                ##      ( pars$m$stagemean[yi,typ] -                      
                ##       sum( sapply(prevtyp,function(ptyp){
                ##         abs(devels[ptyp])^-1 *
                ##           ( pars$m$stagemean[yi,ptyp] +
                ##            pars$m$trend[yi]*devels[ptyp] ) /
                ##              sum(abs(devels)^-1)
                ##       }) ) ) -
                ##      ## sum( sapply( c(1:info$nclus$mir)[-iclus],function(omir){
                ##      ##   ( pars$interactions$cluster[yi,iclus] *
                ##      ##    pars$interactions$cluster[yi,omir] +
                ##      ##    intClusPrecInv[[yi]][iclus,omir] ) *
                ##      ##      pars$mir$stagemean[omir,typ]
                ##      ## }) ) - 
                ##      sum( ( pars$interactions$cluster[yi,iclus] *
                ##            pars$interactions$cluster[yi,omirind] +
                ##            intClusPrecInv[[yi]][iclus,omirind] ) *
                ##          pars$mir$stagemean[omirind,typ] ) - 
                ##      sum( sapply(prevtyp,function(ptyp){
                ##        sum( ( pars$interactions$cluster[yi,iclus] *
                ##              pars$interactions$cluster[yi,] +
                ##              intClusPrecInv[[yi]][iclus,] ) *
                ##            (-pars$mir$stagemean[,ptyp]) ) *
                ##              abs(devels[ptyp])^-1 / sum(abs(devels)^-1)
                ##      }) ) )
                ## })

                allmeans.aschild <- wtdAvgInvDev * expect$m$stageprec *
                  ( pars$interactions$cluster[,iclus] *
                   ( pars$m$stagemean[,typ] -                      
                    apply( sapply(prevtyp,function(ptyp){
                      abs(devels[ptyp])^-1 *
                        ( pars$m$stagemean[,ptyp] +
                         pars$m$trend*devels[ptyp] ) /
                           sum(abs(devels)^-1)
                    }),1,sum) ) -
                   ( matrix(pars$interactions$cluster[,iclus],
                            nrow=info$nclus$m,ncol=length(omirind)) *
                    pars$interactions$cluster[,omirind] +
                    icpi.array[,iclus,omirind] ) %*%
                   pars$mir$stagemean[omirind,typ] - 
                   apply( sapply(prevtyp,function(ptyp){
                     ( matrix(pars$interactions$cluster[,iclus],
                              nrow=info$nclus$m,ncol=info$nclus$mir) *
                      pars$interactions$cluster +
                      icpi.array[,iclus,] ) %*%
                        (-pars$mir$stagemean[,ptyp]) *
                          abs(devels[ptyp])^-1 / sum(abs(devels)^-1)
                   }),1,sum) )

                meansum.aschild <- sum(allmeans.aschild)
              }

              posttyp <- info$typeisbefore$m[[typ]]
              if( length(posttyp)==0 ){
                meansum.asparent <- 0
              }else{
                avgPostWts <- sapply( posttyp,function(ptyp){
                  devels <- pars$m$development[[ptyp]]
                  sum( abs(devels)^-2) / sum( abs(devels)^-1 )
                })

                ## ## old, slow
                ## {Rprof()
                ## allmeans.asparent1 <- sapply( posttyp,function(ptyp){
                ##   sapply( 1:info$nclus$m,function(yi){
                ##     allparents <- info$typefollows[[ptyp]]
                ##     parentmirmat <- pars$mir$stagemean[,allparents,drop=FALSE]
                ##     parentmirmat[iclus,typ] <- 0 # this is the value we're estimating
                ##     devels <- pars$m$development[[ptyp]]
                ##     ret <- avgPostWts[ptyp] * expect$m$stageprec[yi] *
                ##       abs(devels)[typ]^-1 / sum(abs(devels)^-1) *
                ##         -( pars$interactions$cluster[yi,iclus] *
                ##           ( pars$m$stagemean[yi,ptyp] -
                ##            sum( abs(devels)^-1 *
                ##                ( pars$m$stagemean[yi,allparents] +
                ##                 pars$m$trend[yi]*devels) ) /
                ##            sum( abs(devels)^-1 ) ) +
                ##           sum( sapply(allparents,function(xps){
                ##             ( pars$interactions$cluster[yi,iclus] *
                ##              pars$interactions$cluster[yi,] +
                ##              intClusPrecInv[[yi]][iclus,] ) *
                ##                ( -pars$mir$stagemean[,ptyp] +
                ##                 parentmirmat[,xps] ) * 
                ##                   abs(devels[xps])^-1 / sum(abs(devels)^-1)
                ##           }) ) )
                ##     return(ret)
                ##   })
                ## })
                ##  Rprof(NULL)}

                allmeans.asparent <- sapply( posttyp,function(ptyp){
                  allparents <- info$typefollows$m[[ptyp]]
                  parentmirmat <- pars$mir$stagemean[,allparents,drop=FALSE]
                  parentmirmat[iclus,typ] <- 0 # this is the value we're estimating
                  devels <- pars$m$development[[ptyp]]
                  ret <- avgPostWts[ptyp] * expect$m$stageprec *
                    abs(devels)[typ]^-1 / sum(abs(devels)^-1) *
                      -( pars$interactions$cluster[,iclus] *
                        ( pars$m$stagemean[,ptyp] -
                         apply( matrix(abs(devels)^-1,nrow=info$nclus$m,
                                       ncol=length(devels),byrow=TRUE) *
                               ( pars$m$stagemean[,allparents,drop=FALSE] +
                                matrix(pars$m$trend,nrow=info$nclus$m,
                                       ncol=length(devels),byrow=FALSE) *
                                matrix(devels,nrow=info$nclus$m,
                                       ncol=length(devels),byrow=TRUE) ),1,sum) /
                         sum( abs(devels)^-1 ) ) +
                        apply( sapply(allparents,function(xps){
                          ( matrix(pars$interactions$cluster[,iclus],
                                   nrow=info$nclus$m,ncol=info$nclus$mir) *
                           pars$interactions$cluster +
                           icpi.array[,iclus,] ) %*%
                             ( -pars$mir$stagemean[,ptyp] +
                              parentmirmat[,xps] ) * 
                                abs(devels[xps])^-1 / sum(abs(devels)^-1)
                        }),1,sum) )
                  return(ret)
                })
                
                meansum.asparent <- sum(allmeans.asparent)
              }
              meansum <- meansum.aschild + meansum.asparent
            }else{
              meansum <- 0
            }
            return(meansum)
          })
          

          ##cat( interpart/interpart.prec, "\n")
          
          ## put the parts together
          
          ## allprecsum <- diag(as.vector(priorpart.prec)) +
          ##   diag(as.vector(postpart.prec)) + diag(as.vector(cluspart.prec)) +
          ##     interpart.prec[[1]]
          
          pars[[side]]$stagemean.prec[iclus,dotypes] <- priorpart.prec +
            postpart.prec + cluspart.prec + interpart.prec
          pars[[side]]$stagemean[iclus,dotypes] <- 
            ( priorpart + postpart + cluspart + interpart ) /
              pars[[side]]$stagemean.prec[iclus,dotypes]
          
        }
    }
  
  return(pars)
}


updateStageprec <- function( xpr, pars, priors, info, side=NULL,
                            forllh=FALSE, skipcalc=FALSE,
                            intClusPrecInv=NULL){
  expect <- getExpectations( pars )
  ## pars <- updateStageprecBTAhyp( xpr, pars, priors, info, side=side )

  ## stages which utilize *stageprec* as a precision
  usestages <- which( sapply(info$typefollows[[side]],function(x)
                             !("prior" %in% x)) )
  usetypes <- info$types[usestages]

  if( length(usestages)>0 ){
    
    mireffect.list <- getMiReffect(pars,info,side=side)

    if( skipcalc ){ ## this skips the calculation if side="mir",
      mireffectVAR.list <- lapply(info$types,function(typ)
                                  rep(0,info$nclus[[side]]))
    }else{
      mireffectVAR.list <- getMiReffectVAR(pars,info,side=side,
                                           intClusPrecInv=intClusPrecInv)
    }

    ## #######################
    ## the expected value for each stage
    priormeans <- sapply(usetypes,function(typ){
      prevtyp <- info$typefollows[[side]][[typ]]
      devels <- pars[[side]]$development[[typ]]
      ## the weighted average of preceding stages plus the mir targeting
      ## effect by inverse development
      prevmeans <- sapply( 1:info$nclus[[side]], function(xi){
        sum( (pars[[side]]$stagemean[xi,unlist(prevtyp)] +
              mireffect.list[[typ]][xi,] +
              pars[[side]]$trend[xi]*devels ) *
            abs(devels)^-1 ) / sum(abs(devels)^-1)
      })
      ret <- prevmeans
      return(ret)
    })

    ## #######################
    ## calculate the variances associated with above *prevmeans*
    priormean.var <- sapply(usetypes,function(typ){
      prevtyp <- info$typefollows[[side]][[typ]]
      devels <- pars[[side]]$development[[typ]]
      ## the weighted average of variances by inverse development
      prevmean.var <- sapply( 1:info$nclus[[side]], function(xi){
        mireffectVAR.list[[typ]][xi] +
          sum( ( pars[[side]]$stagemean.prec[xi,unlist(prevtyp)]^-1 +
                pars[[side]]$trend.prec[xi]^-1 * devels^2 ) *
              abs(devels)^-2 ) / sum(abs(devels)^-1)^2
      })
      ret <- prevmean.var
      return(ret)
    })

    ## the weighted average of developments of preceding stages
    allWTs <- sapply(usetypes,function(typ){
      devels <- pars[[side]]$development[[typ]]
      avgWTS <- sapply( 1:info$nclus[[side]], function(xi){
        sum( abs(devels)^-2 ) / sum(abs(devels)^-1)
      })
      return(avgWTS)
    })

    alfpart <- rep( length(info$types[usestages]), info$nclus[[side]] )
    btapart <- apply( allWTs * ( (priormeans - pars[[side]]$stagemean[,usestages])^2 +
                                priormean.var + 1/pars[[side]]$stagemean.prec[,usestages] ),
                     1, sum )
  }else{
    alfpart <- rep( 0, info$nclus[[side]] )
    btapart <- rep( 0, info$nclus[[side]] )
  }
  
  pars[[side]]$stageprec$alf <- priors[[side]]$stageprec$alf + 0.5*alfpart
  pars[[side]]$stageprec$bta <- expect[[side]]$stageprecBTAhyp + 0.5*btapart

  if( forllh ){
    ret <- list(psmeans=priormeans,
                psmean.vars=priormean.var,
                allWTs=allWTs)
  }else{
    ret <- pars
  }
  return(ret)
}


updateStageprecBTAhyp <- function( xpr, pars, priors, info, side=NULL ){
  expect <- getExpectations( pars )

  pars[[side]]$stageprecBTAhyp$alf <- priors[[side]]$stageprecBTAhyp$alf +
    priors[[side]]$stageprec$alf * length(expect[[side]]$stageprec)
  pars[[side]]$stageprecBTAhyp$bta <- priors[[side]]$stageprecBTAhyp$bta +
    sum( expect[[side]]$stageprec )

  return(pars)
}


getMembership <- function( xpr,pars,priors,info,
                          side=NULL,shuffle=TRUE,
                          intClusPrecInv=NULL){

  if( any(priors[[side]]$clusmship==0) ){
    pars[[side]]$clusmship <- priors[[side]]$clusmship
  }else if( info$nclus[[side]]==length(info$genes[[side]]) ){
    ## do nothing
  }else{
    
    expect <- getExpectations( pars )
    mship <- matrix(0, length(info$genes[[side]]), info$nclus[[side]] )

    ## #########################
    ## Bayesian clustering
    temp <- t(sapply(info$genes[[side]],function(gn){
      slidemeans <- pars[[side]]$slidemean[gn,]
      slidemeanprecs <- pars[[side]]$slidemean.prec[gn,]
      temprow <- sapply(1:info$nclus[[side]],function(xi){
        expprec <- expect[[side]]$clusprec[xi,]
        stagemeans <- pars[[side]]$stagemean[xi,]
        stagemeanprecs <- pars[[side]]$stagemean.prec[xi,]
        rowitem <- (-0.5*length(stagemeans))*log(2*pi) +
          0.5* sum(gamExpLogx(pars[[side]]$clusprec$alf[xi,], 
                              pars[[side]]$clusprec$bta[xi,])) - 
                                0.5 * ( sum(slidemeans^2*expprec) + sum(expprec/slidemeanprecs) +
                                       sum(stagemeans^2*expprec) + sum(expprec/stagemeanprecs) -
                                       2*sum(slidemeans*expprec*stagemeans) )
        return(rowitem)
      })
      return(temprow)
    }))

    logllh <- log(priors[[side]]$clusmship) + temp

    mship <- t(apply(logllh,1,function(llhrow){
      llhrow[is.na(llhrow)] <- -Inf
      grow <- exp(llhrow-max(llhrow,na.rm=TRUE))
      ## cat(llhrow , "\n")
      ## cat(grow , "\n")
      if( any(is.na(grow)) ){
        ## ## ind <- which(grow==max(grow))
        ## ## retrow <- rep(0,length(grow))
        ## ## retrow[ind] <- 1/length(ind)
        ## grow <- rgamma(length(grow),1)
        grow[is.na(grow)] <- 0  ##/sum(grow)
        retrow <- grow
      }else if( sum(grow)==0 ){
        retrow <- grow
        retrow[] <- 0
      }else{
        retrow <- grow/sum(grow)
      }
      ##retrow[retrow<1e-10] <- 0
      return(retrow)
    }))

    ## pars[[side]]$clusmship <- mship
    
    if( shuffle ){

      ## if a cluster is empty
      clustmemb <- apply(mship>0.1,2,sum)
      emptyclusters <- which( clustmemb==0  )
      if( length(emptyclusters) != 0 ){
        numemptyclusters <- length(emptyclusters)
        cat( " --There are ",numemptyclusters," empty clusters.\n" )
      }

      for( xi in emptyclusters ){
        clustmemb <- apply(mship>0.1,2,sum)
        bigclusts <- which(clustmemb>mean(clustmemb))
        bigclustermemb <- which(apply(mship[,bigclusts,drop=FALSE],1,max)>0.1) #find(max(mship[biggestclust,],[],1)>0.1)
        maxlogllh <- apply(logllh,1,max) #max(logllh)
        ##loneliestgene <- find( maxlogllh == min(maxlogllh(clustmemb==max(clustmemb))) )
        loneliestgene <- which(maxlogllh==min(maxlogllh[bigclustermemb]))

        mship[loneliestgene,] <- 0
        mship[loneliestgene,xi] <- 1
        pars[[side]]$stagemean[xi,] <- pars[[side]]$slidemean[loneliestgene,]
        pars[[side]]$stagemean.prec[xi,] <- max(pars[[side]]$stagemean.prec[-xi,])
        pars[[side]]$clusprec$alf[xi,] <- max(pars[[side]]$clusprec$alf[-xi,])
        pars[[side]]$clusprec$bta[xi,] <- min(pars[[side]]$clusprec$bta[-xi,])
        
        maxlogllh[loneliestgene] <- max(maxlogllh)+1
      }

      ## if two clusters converge
      clusdists <- dist( pars[[side]]$stagemean )
      if( any(clusdists < 1e-5*median(clusdists)) ){
        cat( " --There are some very similar clusters.\n" )

        ind <- which( as.matrix(clusdists) < 0.05*median(clusdists) &
                     as.matrix(clusdists) != 0,
                     arr.ind=TRUE)
        sameclusts <- ind[1,]
        reloclust <- sameclusts[1]

        maxlogllh <- apply(logllh,1,max) #max(logllh)
        loneliestgene <- which.min(maxlogllh)

        ## put the "loneliest" gene into one of the clusters
        mship[loneliestgene,] <- 0
        mship[loneliestgene,reloclust] <- 1
        pars[[side]]$stagemean[reloclust,] <- pars[[side]]$slidemean[loneliestgene,]
        pars[[side]]$stagemean.prec[reloclust,] <- max(pars[[side]]$stagemean.prec[-reloclust,])
        pars[[side]]$clusprec$alf[reloclust,] <- max(pars[[side]]$clusprec$alf[-reloclust,])
        pars[[side]]$clusprec$bta[reloclust,] <- min(pars[[side]]$clusprec$bta[-reloclust,])
        
        maxlogllh[loneliestgene] <- max(maxlogllh)+1
      }
      
      pars[[side]]$clusmship <- mship
      pars <- dataSideUpdate( xpr, pars, priors,info,
                             side=side,shuffle=FALSE,
                             partial=TRUE,
                             intClusPrecInv=intClusPrecInv )
    }

  }
  
  return(pars)
}


updateClusprecs <- function( xpr, pars, priors, info, side=NULL ){
  expect <- getExpectations( pars )

  nummembs <- apply(pars[[side]]$clusmship,2,sum)
  pars[[side]]$clusprec$alf <- priors[[side]]$clusprec$alf +
    0.5 * sapply(info$types,function(typ) nummembs)
  pars[[side]]$clusprec$bta <- expect[[side]]$clusprecBTAhyp +
    0.5 * sapply(info$types,function(typ){
      sapply(1:info$nclus[[side]],function(xi){
        sum( (pars[[side]]$slidemean[,typ]^2 + 1/pars[[side]]$slidemean.prec[,typ] +
              pars[[side]]$stagemean[xi,typ]^2 + 1/pars[[side]]$stagemean.prec[xi,typ] -
              2*pars[[side]]$slidemean[,typ]*pars[[side]]$stagemean[xi,typ]) *
            pars[[side]]$clusmship[,xi] )
      })
    })
  
  return(pars)
}

updateClusprecBTAhyp <- function( xpr, pars, priors, info, side=NULL ){
  expect <- getExpectations( pars )

  pars[[side]]$clusprecBTAhyp$alf <- priors[[side]]$clusprecBTAhyp$alf +
    priors[[side]]$clusprec$alf * prod(dim(expect[[side]]$clusprec))
  pars[[side]]$clusprecBTAhyp$bta <- priors[[side]]$clusprecBTAhyp$bta +
    sum( expect[[side]]$clusprec )

  return(pars)
}


updateTrend <- function( xpr, pars, priors, info, side=NULL ){
  expect <- getExpectations( pars )

  usestages <- which( sapply(info$typefollows[[side]],function(x)
                             !("prior" %in% x)) )
  usetypes <- info$types[usestages]

  if( length(usetypes)>0 ){

    allprecs <- sapply(1:info$nclus[[side]],function(xi){
      sum( sapply(usetypes,function(typ){
        devels <- pars[[side]]$development[[typ]]
        wtdAvgInvDev <- sum(abs(devels)^-2) / sum( abs(devels)^-1 )
        ret <- expect[[side]]$stageprec[xi] * wtdAvgInvDev *
          sum(sign(devels))^2 / sum(abs(devels)^-1)^2 
        return(ret)
      }) )
    })

    mireffect.list <- getMiReffect(pars,info,side=side)

    meansums <- sapply(1:info$nclus[[side]],function(xi){
      sum( sapply(usetypes,function(typ){
        devels <- pars[[side]]$development[[typ]]
        wtdAvgInvDev <- sum(abs(devels)^-2) / sum( abs(devels)^-1 )
        ret <- expect[[side]]$stageprec[xi] * wtdAvgInvDev *
          ( pars[[side]]$stagemean[xi,typ] -
           sum( abs(devels)^-1 *
               ( pars[[side]]$stagemean[xi,info$typefollows[[side]][[typ]]] +
                mireffect.list[[typ]][xi,] ) ) / sum( abs(devels)^-1 )
           ) * sum(sign(devels)) / sum( abs(devels)^-1 )
        return(ret)
      }) )
    })
    
  }else{
    allprecs <- rep(1e10,info$nclus[[side]])
    meansums <- rep(0,info$nclus[[side]])
  }
  
  pars[[side]]$trend.prec <- priors[[side]]$trend$prec + allprecs
  pars[[side]]$trend <- ( priors[[side]]$trend$mean + meansums ) /
    pars[[side]]$trend.prec
  
  return(pars)
}


updateInteractionsCluster <- function( xpr, pars, priors, info ){
  expect <- getExpectations( pars )

  usestages <- which( sapply(info$typefollows$m,function(x)
                             !("prior" %in% x)) )
  usetypes <- info$types[usestages]
  
  allprecs <- lapply(1:info$nclus$m,function(xi){
    lapply(usetypes,function(typ){
      devels <- pars$m$development[[typ]]
      wtdAvgInvDev <- sum(abs(devels)^-2) /
        sum( abs(devels)^-1 )
      prevtyp <- info$typefollows$m[[typ]]
      mirchgavg <- sapply( 1:info$nclus$mir, function(yi){
        sum( abs(devels)^-1 *
            ( pars$mir$stagemean[yi,typ] -
             pars$mir$stagemean[yi,prevtyp] ) ) /
               sum( abs(devels)^-1 )
      })
      mirchgavg.var <- sapply( 1:info$nclus$mir, function(yi){
        sum( abs(devels)^-2 *
            ( pars$mir$stagemean.prec[yi,typ]^-1 +
             pars$mir$stagemean.prec[yi,prevtyp]^-1 ) ) /
               sum( abs(devels)^-1 )^2
      })
      ##mirchgavg.var[] <- 0
      ret <- wtdAvgInvDev * expect$m$stageprec[xi] *
        ( as.matrix(mirchgavg) %*% t(mirchgavg) + diag(mirchgavg.var) )
      
      return(ret)
    })
  })
  precsums <- lapply(allprecs,function(x){
    tempsum <- rep(0,info$nclus$mir)
    for( i in 1:length(x) ){
      ##for( j in 1:length(x[[i]]) ){
      tempsum <- tempsum + x[[i]] #[[j]]
      ##}
    }
    return(tempsum)
  })

  
  allmeans <- lapply(1:info$nclus$m,function(xi){
    lapply(usetypes,function(typ){
      devels <- pars$m$development[[typ]]
      wtdAvgInvDev <- sum(abs(devels)^-2) / sum( abs(devels)^-1 )
      prevtyp <- info$typefollows$m[[typ]]
      mirchgavg <- sapply( 1:info$nclus$mir, function(yi){
        sum( abs(devels)^-1 *
            ( pars$mir$stagemean[yi,typ] -
             pars$mir$stagemean[yi,prevtyp] ) ) /
               sum( abs(devels)^-1 )
      })
      ret <- wtdAvgInvDev * expect$m$stageprec[xi] *
        mirchgavg * ( pars$m$stagemean[xi,typ]  -
                     sum( abs(devels)^-1 *
                         ( pars$m$stagemean[xi,prevtyp] +
                          pars$m$trend[xi] * devels ) /
                         sum(abs(devels)^-1) ) )
      return(ret)
    })
  })
  meansums <- lapply(allmeans,function(x){
    tempsum <- rep(0,info$nclus$mir)
    for( i in 1:length(x) ){
      tempsum <- tempsum + as.vector(x[[i]])
    }
    return(tempsum)
  })

  
  ## the sum of individual interactions, the prior
  indivsum <- t(pars$m$clusmship) %*%
    pars$interactions$individual %*%
      pars$mir$clusmship

  indivsum2 <- sapply( 1:info$nclus$mir,function(mirclus){
    sapply( 1:info$nclus$m,function(mclus){
      sum( pars$interactions$individual *
          ( pars$m$clusmship[,mclus,drop=FALSE] %*%
           t(pars$mir$clusmship[,mirclus,drop=FALSE]) ) )
    })
  })
  
  ## indiv.members <- sapply( 1:info$nclus$mir,function(mirclus){
  ##   sapply( 1:info$nclus$m,function(mclus){
  ##     sum( pars$m$clusmship[,mclus,drop=FALSE] %*%
  ##         t(pars$mir$clusmship[,mirclus,drop=FALSE]) )
  ##   })
  ## })
  ## 
  ## postpart <- indivsum * expect$interactions$prec.common
  ## postpart.prec <- expect$interactions$prec.common * indiv.members

  pars$interactions$cluster.prec <- lapply(1:info$nclus$m,function(xi){
    diag(expect$interactions$prec.common*rep(1,info$nclus$mir)) + precsums[[xi]]
  })
  pars$interactions$cluster <- t(sapply(1:info$nclus$m,function(xi){
    ginv( pars$interactions$cluster.prec[[xi]] ) %*% 
      ( expect$interactions$prec.common * indivsum[xi,] + meansums[[xi]] )
  }) )
  
  return(pars)
}


getParmat <- function(mname,mirname,possmirs,info){
  ind <- which(possmirs[[mname]][,"miRNA.main"]==mirname)
  if( length(ind)>0 ){
    parmat <- possmirs[[mname]][ind,info$interactionfactors,
                                drop=FALSE]
    parmat[parmat=="NULL"] <- "0"
    ##parmat <- apply(parmat,2,as.numeric)
    parmat <- apply(parmat,2,function(x) sum(as.numeric(x)))
    if(is.null(dim(parmat))) parmat <- t(parmat)
  }else{
    ##parmat <- pars$interactions$intercoeffs
    ##parmat[] <- 0
    parmat <- matrix(0,1,length(info$interactionfactors))
    colnames(parmat) <- info$interactionfactors
  }
  parmat[,"constant"] <- 1 #pmin(parmat[,"constant"],1)

  return(parmat)
}

buildParmatList <- function(possmirs,info){
  parmat.list <- lapply(info$genes$mir,function(mirgene){
    mirname <- info$possmir.genes$mir[mirgene]
    parmat <- sapply(info$genes$m,function(mgene){
      mname <- info$possmir.genes$m[mgene]
      ret <- getParmat(mname,mirname,possmirs,info)
    } )
    if( is.null(dim(parmat)) ){
      parmat <- as.matrix(parmat)
    }else{
      parmat <- t(parmat)
    }
    colnames(parmat) <- info$interactionfactors
    rownames(parmat) <- info$genes$m
    return(parmat)
  })
  names(parmat.list) <- info$genes$mir
  return(parmat.list)
}


updateInteractionsIndividual <- function(xpr, pars, priors, info,
                                         parmat.list){
  expect <- getExpectations( pars )

  t.m.mship <- t(pars$m$clusmship)


  for( mirgene in info$genes$mir ){

    mirname <- info$possmir.genes$mir[mirgene]
    priorpart.prec <- info$prediction.weight * expect$interactions$prec.common
    parmat <- parmat.list[[mirgene]]

    priorpart.all <- info$prediction.weight *
      expect$interactions$prec.common *
        t(pars$interactions$intercoeffs) %*% t(parmat)

    for( mgene in info$genes$m ){

      priorpart <- priorpart.all[,mgene,drop=FALSE]
      
      ## bicluspair <- t(pars$m$clusmship[mgene,,drop=FALSE]) %*%
      ##   pars$mir$clusmship[mirgene,,drop=FALSE]
      bicluspair <- t.m.mship[,mgene,drop=FALSE] %*%
        pars$mir$clusmship[mirgene,]

      postpart.prec <- expect$interactions$prec.common * sum( bicluspair )
      
      
      useindivinters <- pars$interactions$individual
      useindivinters[mgene,mirgene] <- 0
      ## other.contribs <- sapply( 1:info$nclus$mir,function(mirclus){
      ##   sapply( 1:info$nclus$m,function(mclus){
      ##     sum( useindivinters * 
      ##         ( pars$m$clusmship[,mclus,drop=FALSE] %*%
      ##          t(pars$mir$clusmship[,mirclus,drop=FALSE]) ) )
      ##   })
      ## })
      ## other.contribs <- sapply( 1:info$nclus$mir,function(mirclus){
      ##   sapply( 1:info$nclus$m,function(mclus){
      ##     t.m.mship[mclus,,drop=FALSE] %*%
      ##       ##t(pars$m$clusmship[,mclus,drop=FALSE]) %*%
      ##       useindivinters %*% 
      ##         pars$mir$clusmship[,mirclus,drop=FALSE]
      ##   })
      ## })

      other.contribs <- t.m.mship %*%
        useindivinters %*% pars$mir$clusmship

      
      postpart <- expect$interactions$prec.common *
        sum( bicluspair * ( pars$interactions$cluster - other.contribs ) )
      

      pars$interactions$individual.prec[mgene,mirgene] <-
        priorpart.prec + postpart.prec
      pars$interactions$individual[mgene,mirgene] <- (priorpart + postpart) /
        (priorpart.prec + postpart.prec)
      
      ## in case of weird data from *possmirs*
      pars$interactions$individual[is.na(pars$interactions$individual)] <- 0
    }
    
  }
  
  return(pars)
}



updateInteractionsCoefficients <- function(xpr, pars, priors, info, parmat.list){
  expect <- getExpectations( pars )

  allprecs <- lapply(info$genes$mir,function(mirgene){
    lapply(info$genes$m,function(mgene){
      parmat <- parmat.list[[mirgene]][mgene,,drop=FALSE]
      prec <- info$prediction.weight *
        expect$interactions$prec.common *
          t(parmat) %*% parmat
      return(prec)
    })
  })
  
  allmeans <- lapply(info$genes$mir,function(mirgene){
    lapply(info$genes$m,function(mgene){
      parmat <- parmat.list[[mirgene]][mgene,,drop=FALSE]

      ret <- info$prediction.weight *
        expect$interactions$prec.common *
          pars$interactions$individual[mgene,mirgene] *
            parmat
      
      return(ret)
    })
  })

  precsum <- matrix(0,length(info$interactionfactors),
                    length(info$interactionfactors))
  for( i in 1:length(allprecs) ){
    for( j in 1:length(allprecs[[i]]) ){
      precsum <- precsum + allprecs[[i]][[j]]
    }
  }

  meansum <- rep(0,length(info$interactionfactors))
  for( i in 1:length(allmeans) ){
    for( j in 1:length(allmeans[[i]]) ){
      meansum <- meansum + allmeans[[i]][[j]]
    }
  }

  pars$interactions$intercoeffs.prec <- priors$interactions$intercoeffs$prec +
    precsum
  pars$interactions$intercoeffs <-
    as.vector(ginv(pars$interactions$intercoeffs.prec) %*%
              t( priors$interactions$intercoeffs$mean + meansum )) #qr.solve(precsum,meansum)
  names(pars$interactions$intercoeffs) <- names(priors$interactions$intercoeffs$mean)
  
  return(pars)
}


updateInteractionPrecCommon <- function(xpr, pars, priors, info, parmat.list){
  expect <- getExpectations( pars )
  
  pred2ind.priormeans <- sapply(info$genes$mir,function(mirgene){
    sapply(info$genes$m,function(mgene){
      parmat <- parmat.list[[mirgene]][mgene,,drop=FALSE]
      ret <-  t(pars$interactions$intercoeffs) %*% t(parmat)
      return(ret)
    })
  })

  tempPrecInv <- ginv(pars$interactions$intercoeffs.prec)
  pred2ind.priormeanVARs <- sapply(info$genes$mir,function(mirgene){
    sapply(info$genes$m,function(mgene){
      parmat <- parmat.list[[mirgene]][mgene,,drop=FALSE]
      ret <- sum(diag( t(parmat) %*% parmat %*%
                      tempPrecInv ))
      return(ret)
    })
  })

  clustVARs <- t(sapply(pars$interactions$cluster.prec,function(x){
    diag(ginv(x))
  }))

  ind2clust.priormeans <- t(pars$m$clusmship) %*%
    pars$interactions$individual %*% pars$mir$clusmship

  ind2clust.priormeanVARs <- t(pars$m$clusmship) %*%
    (1/pars$interactions$individual.prec) %*% pars$mir$clusmship

  
  ## temp <- sapply(1:info$nclus$mir,function(mirclus){
  ##   sapply(1:info$nclus$m,function(mclus){
  ##     membmat <- pars$m$clusmship[,mclus,drop=FALSE] %*%
  ##       t(pars$mir$clusmship[,mirclus,drop=FALSE])
  ##     membintvec <- as.vector(pars$interactions$individual *
  ##                             membmat)
  ##     pairwise.prods <- as.matrix(membintvec) %*% t(membintvec)
  ##     diag(pairwise.prods) <- 0
  ##     ret <- pars$interactions$cluster[mclus,mirclus]^2 +
  ##       clustVARs[mclus,mirclus] +
  ##         sum( (pars$interactions$individual^2 +
  ##               1/pars$interactions$individual.prec) * membmat ) +
  ##                 sum( pairwise.prods ) -
  ##                   2 * sum( pars$interactions$individual * membmat *
  ##                           pars$interactions$cluster[mclus,mirclus] )              
  ##     return(ret)
  ##   })
  ## })

  
  pars$interactions$prec.common$alf <- priors$interactions$prec.common$alf +
    0.5 * prod(dim(pars$interactions$individual)) +
      0.5 * prod(dim(pars$interactions$cluster))
  
  pars$interactions$prec.common$bta <- priors$interactions$prec.common$bta +
    0.5 * info$prediction.weight *
      sum( (pars$interactions$individual - pred2ind.priormeans)^2 +
          1/pars$interactions$individual.prec + pred2ind.priormeanVARs ) +
            0.5 * sum( (pars$interactions$cluster - ind2clust.priormeans)^2 +
                      clustVARs + ind2clust.priormeanVARs )
  
  return(pars)
}



updateAllInteractions <- function(xpr, pars, priors, info, parmat.list){

  ## for( i in 1:10 ){

  ## llh <- list()
  ## llh[[ length(llh)+1 ]] <- calcllh(priors,pars,info,possmirs)
    
  pars <- updateInteractionPrecCommon(xpr, pars, priors, info, parmat.list)
  pars <- updateInteractionsCluster( xpr, pars, priors, info )
  pars <- updateStageprec( xpr, pars, priors, info, side="m" )
  pars <- updateInteractionsIndividual(xpr, pars, priors, info, parmat.list)
  pars <- updateInteractionsCoefficients(xpr, pars, priors, info, parmat.list)
    
  ## cat(diff(sapply(llh,function(x) x[[1]]))) ##,"\n")
  ## indivsum <- t(pars$m$clusmship) %*%
  ##   pars$interactions$individual %*%
  ##     pars$mir$clusmship
  ## cat(" -- ", sum( (pars$interactions$cluster - indivsum)^2 ),
  ##     pars$interactions$prec.common$bta,"\n")
  ## }
  
  return(pars)
}


initializeInteractions <- function(xpr, pars, priors, info,
                                   possmirs, parmat.list){

  pars$interactions$prec.common$alf[] <- 1e-1
  pars$interactions$prec.common$bta[] <- 1
  pars$m$stageprec$alf[] <- 1e5
  pars$m$stageprec$bta[] <- 1e-5

  ## pars$interactions$individual <- sapply(info$genes$mir,function(mirgene){
  ##   sapply(info$genes$m,function(mgene){
  ##     mirname <- info$possmir.genes$mir[mirgene]
  ##     mname <- info$possmir.gene$m[mgene]
  ##     ind <- which(possmirs[[mname]][,"miRNA.main"]==mirname)
  ##     return(-length(ind)/5)
  ##   })
  ## })

  ## savedstuff <- list(xpr=xpr,pars=pars,priors=priors,
  ##                    info=info,
  ##                    possmirs=possmirs)
  ## save(savedstuff,file="savedstuff.Rdata")

  cat("Initial updates to interaction variables...\n")
  for( i in 1:3 ){
    pars <- updateInteractionsCluster( xpr, pars, priors, info )
    ##pars <- updateInteractionPrecCluster(xpr, pars, priors, info, parmat.list)
    pars <- updateInteractionPrecCommon(xpr, pars, priors, info, parmat.list)
    pars <- updateStageprec( xpr, pars, priors, info, side="m" )
  }

  ## savedstuff <- list(xpr=xpr,pars=pars,priors=priors,
  ##                    info=info,
  ##                    possmirs=possmirs)
  ## save(savedstuff,file="savedstuff.Rdata")

  for( i in 1:3 ){
    pars <- updateInteractionsIndividual(xpr, pars, priors, info, parmat.list)
    pars <- updateInteractionsCoefficients(xpr, pars, priors, info, parmat.list)
    ##pars <- updateInteractionPrecIndividual(xpr, pars, priors, info, parmat.list)
    pars <- updateInteractionPrecCommon(xpr, pars, priors, info, parmat.list)
  }
  
  return(pars)
}


dataSideUpdate <- function( xpr, pars, priors, info,
                           side=NULL,shuffle=TRUE,
                           partial=FALSE,
                           intClusPrecInv=NULL){
  ## llh <- list()
  ## llh[[ length(llh)+1 ]] <- calcllh(priors,pars,info,possmirs)

  if( !partial ){
    pars <- getMembership( xpr,pars,priors,info,side=side,shuffle=shuffle,
                          intClusPrecInv=intClusPrecInv)
  }

  pars <- updateStages( xpr, pars, priors, info, side=side,
                       intClusPrecInv=intClusPrecInv)
  pars <- updateStageprec( xpr, pars, priors, info, side=side,
                          intClusPrecInv=intClusPrecInv)
  pars <- updateWithinSlide( xpr, pars, priors, info, side=side )
  pars <- updateAcrossSlides( xpr, pars, priors, info, side=side )
  pars <- updateClusprecs( xpr, pars, priors, info, side=side)
  pars <- updateClusprecBTAhyp( xpr, pars, priors, info, side=side )
  pars <- updateStageprecBTAhyp( xpr, pars, priors, info, side=side )
  ## pars <- updateTrend( xpr, pars, priors, info, side=side )

  ## cat(diff(sapply(llh,function(x) x[[1]])),"\n")

  return( pars )
}

initializeDataSide <- function(xpr,pars,priors,info,
                               initialiters=10){

  intClusPrecInv <- lapply(pars$interactions$cluster.prec,ginv)

  for( side in c("mir","m") ){
    ##for( side in c("mir") ){
    for( i in 1:initialiters ){
      pars <- updateWithinSlide( xpr, pars, priors, info, side=side  )
    }
    for( i in 1:initialiters ){
      pars <- updateAcrossSlides( xpr, pars, priors, info, side=side  )
      pars <- updateWithinSlide( xpr, pars, priors, info, side=side  ) 
    }
    pars[[side]]$clusprec$alf[,] <- 10
    pars[[side]]$clusprec$bta[,] <- 1
    cat(side, " initial variable updates... \n")
    
    for( i in 1:round(initialiters/2) ){
      ## pars <- getMembership( xpr,pars,priors,info,side=side,shuffle=TRUE,
      ##                       intClusPrecInv=intClusPrecInv)
      pars <- updateStages( xpr, pars, priors, info, side=side,
                           intClusPrecInv=intClusPrecInv )
    }
    
    cat(side, " expanded variable updates... \n")
    for( i in 1:round(initialiters) ){
      pars <- updateAcrossSlides( xpr, pars, priors, info, side=side  )
      pars <- updateWithinSlide( xpr, pars, priors, info, side=side  )
      pars <- updateClusprecs( xpr, pars, priors, info, side=side)
      ## pars <- getMembership( xpr,pars,priors,info,side=side,shuffle=TRUE,
      ##                       intClusPrecInv=intClusPrecInv)
      pars <- updateStages( xpr, pars, priors, info, side=side,
                           intClusPrecInv=intClusPrecInv)
    }
    
    pars[[side]]$stageprec$alf[] <- 5 #mean(pars[[side]]$clusprec$alf)
    pars[[side]]$stageprec$bta[] <- 1 #2*mean(pars[[side]]$clusprec$bta)
    for( i in 1:round(initialiters/2) ){
      pars <- dataSideUpdate( xpr, pars, priors,info,
                             side=side,shuffle=TRUE,
                             partial=TRUE,
                             intClusPrecInv=intClusPrecInv )
    }
  } ## end loop over *side*
  return(pars)
}


fitModel <- function(xpr,pars,priors,info,possmirs,
                     initialiters=5,
                     dataiters=2,
                     smalldataiters=5,
                     interactioniters=5,
                     bigiters=1,
                     doclustering=FALSE){

  for( side in c("mir","m")){
    if( !doclustering[[side]] ){
      info$nclus[[side]] <- length(info$genes[[side]])
    }
  }
  
  priors <- initpriors(info,possmirs)
  pars <- initpars(info,possmirs)

  for( side in c("mir","m") ){
    if( !doclustering[[side]] ){
      priors[[side]]$clusmship <- diag(rep(1,info$nclus[[side]]))
      rownames(priors[[side]]$clusmship) <- info$genes[[side]]
    }
  }

  parmat.list <- buildParmatList(possmirs,info)

  
  cat("\n")
  timestamp()
  cat("\nBegin model fit with ",info$nclus$mir," mir and ",
      info$nclus$m," m clusters...\n")

  timestamp()
  pars <- initializeDataSide(xpr,pars,priors,info,
                             initialiters=initialiters) 

  timestamp()
  pars <- initializeInteractions(xpr, pars, priors, info,
                                 possmirs, parmat.list)
  

  ##for( bi in 0:(bigiters-1)){
  bi <- 0
  while( bi < bigiters ){
    savedstuff <- list(xpr=xpr,pars=pars,priors=priors,
                       info=info,
                       possmirs=possmirs)
    save(savedstuff,file="savedstuff.Rdata")

    timestamp()
    cat("Big iteration:",bi,"of",bigiters-1,"\n")
    llh0 <- calcllh(priors,pars,info,possmirs)


    ## if( bi == 5 ){
    ##   llh0hyp0 <- calcllh(priors,pars,info,possmirs)
    ##   cat("Opening up some priors to make them less informative before\n ",
    ##       "  optimizing them in a few more big iterations.\n")
    ##   for( side in c("mir","m") ){
    ##     priors[[side]]$clusprecBTAhyp <- list(alf=1e-5,bta=1e-5)
    ##     priors[[side]]$stageprecBTAhyp <- list(alf=1e-5,bta=1e-5)
    ##   }
    ##   llh0hyp1 <- calcllh(priors,pars,info,possmirs)
    ##   cat( "Prior-opening likelihood change (+ or -):",
    ##       llh0hyp1[[1]]-llh0hyp0[[1]],"\n")
    ## }

    if(  bi > 10 ){
      llh0hyp0 <- calcllh(priors,pars,info,possmirs)
      cat("Updating hyper-priors.\n")
      for( i in 1:3 ){
        priors <- updateHyperpriors( xpr, pars, priors, info )
      }
      llh0hyp1 <- calcllh(priors,pars,info,possmirs)
      cat( "Hyper-prior update likelihood improvement:",
          llh0hyp1[[1]]-llh0hyp0[[1]],"\n")
    }
    
    for( di in 0:(dataiters-1) ){  
      timestamp()
      cat("Data iteration:",di,"of",dataiters-1,"\n")
      llh1 <- calcllh(priors,pars,info,possmirs)

      ## pre-calculation
      intClusPrecInv <- lapply(pars$interactions$cluster.prec,ginv)

      for( shuffle in c(TRUE,FALSE) ){
        for( i in 0:(smalldataiters-1) ){
          for( side in c("mir","m") ){

            pars <- dataSideUpdate( xpr, pars, priors, info,
                                   side=side,shuffle=shuffle,
                                   partial=FALSE,
                                   intClusPrecInv=intClusPrecInv)
            pars <- dataSideUpdate( xpr, pars, priors, info,
                                   side=side,shuffle=FALSE,
                                   partial=FALSE,
                                   intClusPrecInv=intClusPrecInv)
          } ## end loop over *side*
        } ## end loop over *smalldataiters*
      } ## end loop over *shuffle*

      if( di != dataiters-1 ){ # if it's not the last data iteration
        ## update fixed parameters
        timestamp()
        cat("updating fixed parameters\n")
        for( side in c("mir","m") ){
          ## update *development
          cat("--development: ")
          for( typ in names(pars[[side]]$development) ){
            cat(typ)
            for( ptyp in names(pars[[side]]$development[[typ]]) ){
              pars <- optimizeDevelopments( typ, ptyp,
                                           priors, pars, info,
                                           side=side,
                                           intClusPrecInv=intClusPrecInv)
              cat(".")
            }
          }
          cat("\n")
        }
      }

      llh2 <- calcllh(priors,pars,info,possmirs)
      cat( "Data iteration likelihood improvement:",
          llh2[[1]]-llh1[[1]],"\n")
    } ## end loop over *dataiters*

    if( interactioniters>0 ){
      llh3 <- calcllh(priors,pars,info,possmirs)
      for( ii in 0:(interactioniters-1)){
        cat("Interaction iteration:",ii,"\n")
        pars <- updateAllInteractions(xpr, pars, priors, info, parmat.list)
      }
      timestamp()
      llh4 <- calcllh(priors,pars,info,possmirs)
      cat( "Interaction iteration likelihood improvement:",
          llh4[[1]]-llh3[[1]],"\n")

      bigiter.llhchg <- llh4[[1]]-llh0[[1]]
      cat( "Big iteration likelihood improvement:",
          bigiter.llhchg,"\n\n")
    }

    llh <- calcllh(priors,pars,info,possmirs)
    ret <- list(llh=llh,pars=pars,
                priors=priors,info=info)
    ##save(ret,file=paste("ret",bi,".Rdata",sep=""))

    bi <- bi+1
    if( abs(bigiter.llhchg) < length(info$genes$mir)/100 &&
       bi > 15 ) bi <- bigiters+1
  } ## end loop over *bigiters*
  
  return(ret)
}

## for fitting with a range of cluster values
trimclusrange <- function(doclusnums,llh.vec,clusrange.total,bestclusnum){
  trimprop <- 0.3
  trimnum <- floor( length(doclusnums)*trimprop )
  n <- length(doclusnums)

  maxind <- which.max(llh.vec)
  if( doclusnums[maxind]==doclusnums[1] ){
    doclusnums <- c(max( max(clusrange.total[1],bestclusnum-1),
                        2*doclusnums[1]-doclusnums[2]),
                    doclusnums[1:(n-trimnum-1)])
  }else if( doclusnums[maxind]==doclusnums[n] ){
    doclusnums <- c(doclusnums[(trimnum+2):n],
                    min( min(clusrange.total[2],bestclusnum+1),
                        2*doclusnums[n]-doclusnums[n-1]))
  }else{  
    for( i in 1:trimnum ){
      newn <- length(doclusnums)
      trimind <- ifelse(llh.vec[1]<llh.vec[newn],1,newn)
      doclusnums <- doclusnums[-trimind]
    }
  }
  newclusrange <- range(doclusnums)
  return(newclusrange)
}


## ###########################################
## non-Bayesian parameters


devoptimfunc <- function(priors,pars,info,side=NULL,
                         intClusPrecInv=NULL){
  llhlist <- list(
                  development=calcllh.development(priors,pars,info,side=side),
                  stagemean=calcllh.stagemean(priors,pars,info,side=side,
                    intClusPrecInv=intClusPrecInv)
                  )
  ret <- sum(unlist(llhlist))
  return(ret)
}

optimizeDevelopments <- function(typ,ptyp,
                                 priors,
                                 pars,
                                 info,
                                 side=NULL,
                                 returninitvals=FALSE,
                                 initvals=NULL,
                                 intClusPrecInv=NULL){

  ## scale.init <- 0.5
  ## val.init <- abs(pars[[side]]$development[[typ]][ptyp])


  out.pos <- optimize( function(x){
    pars[[side]]$development[[typ]][ptyp] <- x
    ret <- devoptimfunc(priors,pars,info,side=side,
                        intClusPrecInv=intClusPrecInv)
    return( ret )
  },c(1e-2,1e1),maximum=TRUE )$maximum

  ## out.pos <- optimize( function(x){
  ##   pars[[side]]$development[[typ]][ptyp] <- x
  ## 
  ##   {Rprof()
  ##   llh <- calcllh(priors,pars,info,possmirs)
  ##    Rprof(NULL)}
  ##    
  ##   ##cat(x,llh[[1]],"\n")
  ##   ret <- llh[[1]]
  ##   return( ret )
  ## },c(1e-2,1e1),maximum=TRUE )$maximum

  ## out.neg <- optimize( function(x){
  ##   pars[[side]]$development[[typ]][ptyp] <- -x
  ##   llh <- calcllh(priors,pars,info,possmirs)
  ##   ##cat(-x,llh[[1]],"\n")
  ##   ret <- llh[[1]]
  ##   return( ret )
  ## },c(1e-1,1e1),maximum=TRUE,tol=0.1 )$maximum

  ## out.pos2 <- metrop( y <- function(xvec){
  ##   vals <- sapply(xvec,function(x){
  ##     if( x <= 0 ){
  ##       ret <- -Inf
  ##     }else{
  ##       pars[[side]]$development[[typ]][ptyp] <- x
  ##       llh <- calcllh(priors,pars,info,possmirs)
  ##       ret <- llh[[1]]
  ##     }
  ##     return( ret )
  ##   })
  ##   return(vals)
  ## },val.init,100,scale=scale.init )$final
  
  ## out.neg <- metrop( y <- function(xvec){
  ##   vals <- sapply(xvec,function(x){
  ##     if( x <= 0 ){
  ##       ret <- -Inf
  ##     }else{
  ##       pars[[side]]$development[[typ]][ptyp] <- -x
  ##       llh <- calcllh(priors,pars,info,possmirs)
  ##       ret <- llh[[1]]
  ##     }
  ##     return( ret )
  ##   })
  ##   return(vals)
  ## },val.init,30,scale=scale.init )$final

  
  ## pars[[side]]$development[[typ]][ptyp] <- out.pos
  ## llh.pos <- calcllh(priors,pars,info,possmirs)[[1]]
  ## pars[[side]]$development[[typ]][ptyp] <- -out.neg
  ## llh.neg <- calcllh(priors,pars,info,possmirs)[[1]]
  ## if( llh.neg > (llh.pos+log(1.5)) ){
  ##   bestval <- -out.neg
  ## }else{
  ##   bestval <- out.pos
  ## }
  ## pars[[side]]$development[[typ]][ptyp] <- bestval

  pars[[side]]$development[[typ]][ptyp] <- out.pos

  return(pars)
}


priorshapederiv <- function(x,parHYPshape,parHYPscale, 
                            parshape,parscale){
  y <- sum( gamExpLogx(parHYPshape,parHYPscale) -
           psigamma(x) +
           gamExpLogx(parshape,parscale) )
  return(y)
}

priorHYPshapederiv <- function(x,priorscale, 
                               parshape,parscale){
  y <- sum( log(priorscale) -
           psigamma(x) +
           gamExpLogx(parshape,parscale) )
  return(y)
}


updateHyperpriors <- function( xpr, pars, priors, info ){
  expect <- getExpectations(pars)


  for( side in c("mir","m") ){

    ## Stagemeans
    stageswithprior <- info$types[ sapply(info$types,function(typ){
      "prior" %in% info$typefollows[[side]][[typ]]
    }) ]

    ## priors[[side]]$stagemean$prec <- length(stageswithprior) *
    ##  info$nclus[[side]] / sum( ( pars[[side]]$stagemean[,stageswithprior] -
    ##                             priors[[side]]$stagemean$mean )^2  +
    ##                           1/pars[[side]]$stagemean.prec[,stageswithprior] )
    
    ## Techprecs
    for( precname in c("spot","slide") ){
      priors[[side]]$techprec[[precname]]$alf <- uniroot( function(x){
        priorHYPshapederiv(x,priors[[side]]$techprec[[precname]]$bta,
                           pars[[side]]$techprec[[precname]]$alf,
                           pars[[side]]$techprec[[precname]]$bta)
      },c(1e-10,1e10) )$root
      priors[[side]]$techprec[[precname]]$bta <-
        priors[[side]]$techprec[[precname]]$alf /
          expect[[side]]$techprec[[precname]]
    }

    ## Clusprec and hyperprior
    priors[[side]]$clusprecBTAhyp$alf <- uniroot( function(x){
      priorHYPshapederiv(x,priors[[side]]$clusprecBTAhyp$bta,
                         pars[[side]]$clusprecBTAhyp$alf,
                         pars[[side]]$clusprecBTAhyp$bta)
    },c(1e-10,1e10) )$root
    priors[[side]]$clusprecBTAhyp$bta <- priors[[side]]$clusprecBTAhyp$alf /
      expect[[side]]$clusprecBTAhyp

    priors[[side]]$clusprec$alf <- uniroot( function(x){
      priorshapederiv(x,pars[[side]]$clusprecBTAhyp$alf,
                      pars[[side]]$clusprecBTAhyp$bta,
                      pars[[side]]$clusprec$alf,
                      pars[[side]]$clusprec$bta)
    },c(1e-10,1e10) )$root

    ## Stageprec and hyperprior
    priors[[side]]$stageprecBTAhyp$alf <- uniroot( function(x){
      priorHYPshapederiv(x,priors[[side]]$stageprecBTAhyp$bta,
                         pars[[side]]$stageprecBTAhyp$alf,
                         pars[[side]]$stageprecBTAhyp$bta)
    },c(1e-10,1e10) )$root
    
    priors[[side]]$stageprecBTAhyp$bta <- priors[[side]]$stageprecBTAhyp$alf /
      expect[[side]]$stageprecBTAhyp

    priors[[side]]$stageprec$alf <- uniroot( function(x){
      priorshapederiv(x,pars[[side]]$stageprecBTAhyp$alf,
                      pars[[side]]$stageprecBTAhyp$bta,
                      pars[[side]]$stageprec$alf,
                      pars[[side]]$stageprec$bta)
    },c(1e-10,1e10) )$root

    ## ## Trend
    ## if( side=="mir" ){
    ##   priors[[side]]$trend$prec <- info$nclus[[side]] /
    ##     sum( ( pars[[side]]$trend -
    ##           priors[[side]]$trend$mean )^2  +
    ##         1/pars[[side]]$trend.prec )
    ## }

  }
  
  ## ###############
  ## INTERACTION

  priors$interactions$prec.common$alf <- uniroot( function(x){
    priorHYPshapederiv(x,priors$interactions$prec.common$bta,
                       pars$interactions$prec.common$alf,
                       pars$interactions$prec.common$bta)
  },c(1e-10,1e10) )$root

  priors$interactions$prec.common$bta <-
    priors$interactions$prec.common$alf /
      mean( expect$interactions$prec.common )

  ncoeffs <- length(priors$interactions$intercoeffs$mean)
  if( ncoeffs>1 ){
    priors$interactions$intercoeffs$prec <-
      diag( rep( ncoeffs /
                sum( ( pars$interactions$intercoeffs -
                      priors$interactions$intercoeffs$mean )^2  +
                    diag(ginv(pars$interactions$intercoeffs.prec)) ),
                ncoeffs) )
  }else{
    priors$interactions$intercoeffs$prec <- ncoeffs /
        sum( ( pars$interactions$intercoeffs -
              priors$interactions$intercoeffs$mean )^2  +
            diag(ginv(pars$interactions$intercoeffs.prec)) )
  }

  
  return(priors)
}

## ##########################################################
## functions co calculate the marginal likelihood lower bound


gamEnt <- function( alf, bta ){
  ent <- alf - log(bta) + lgamma(alf) + (1-alf)*psigamma(alf)
  return(ent)
}

gaussEnt <- function( prec ){
  ent <- 0.5*( 1 + log(2*pi) - log(prec) )
  return(ent)
}

multiGaussEnt <- function( prec ){
  ent <- 0.5*( dim(prec)[1]*log(2*pi*exp(1)) -
              determinant(prec,logarithm=TRUE)$modulus[1] )
  return(ent)
}


## ###########################################
## separate functions for calculating the log likelihood

calcllh.interactions.cluster <- function(priors,pars,info){
  expect <- getExpectations( pars )
  list(
       qlogp=sapply( 1:info$nclus$m,function(mclus){

         ## sqdiffmat2 <- t( sapply( 1:info$nclus$mir,function(mirclus){
         ##   pars$interactions$cluster[mclus,mirclus] -
         ##     sum( pars$m$clusmship[,mclus,drop=FALSE] %*%
         ##         t(pars$mir$clusmship[,mirclus,drop=FALSE]) *
         ##         pars$interactions$individual )
         ## }) )
         sqdiffmat <- pars$interactions$cluster[mclus,,drop=FALSE] -
           t(pars$m$clusmship[,mclus,drop=FALSE]) %*%
             pars$interactions$individual %*%
               pars$mir$clusmship
         
         ## sqdiffmatVAR2 <- ginv(pars$interactions$cluster.prec[[mclus]]) +
         ##   diag(as.vector( sapply( 1:info$nclus$mir,function(mirclus){
         ##     sum( pars$m$clusmship[,mclus,drop=FALSE] %*%
         ##         t(pars$mir$clusmship[,mirclus,drop=FALSE]) *
         ##         1/pars$interactions$individual.prec )
         ##   }) ))

         sqdiffmatVAR <- ginv(pars$interactions$cluster.prec[[mclus]]) +
           diag(as.vector( 
                          t(pars$m$clusmship[,mclus,drop=FALSE]) %*%
                          (1/pars$interactions$individual.prec) %*%
                          pars$mir$clusmship ) )

         precmat <- diag(rep(expect$interactions$prec.common,
                             info$nclus$mir))
         ret <- -0.5 * info$nclus$mir * log(2*pi) +
           0.5 * info$nclus$mir *
             gamExpLogx(pars$interactions$prec.common$alf,
                        pars$interactions$prec.common$bta) -
                          0.5 * ( sqdiffmat %*% precmat %*% t(sqdiffmat) +
                                 sum( diag( precmat * sqdiffmatVAR ) ) )
         return(ret)
       }),
       qlogq=
       -sapply(1:info$nclus$m,function(mclus){
         -multiGaussEnt(pars$interactions$cluster.prec[[mclus]])
       })
       )
}

calcllh.interactions.individual <- function(priors,pars,info){
  expect <- getExpectations( pars )
  intercoeffPrecInv <- ginv(pars$interactions$intercoeffs.prec)
  parmat.list <- buildParmatList(possmirs,info)

  list(
       qlogp=
       sapply(info$genes$mir,function(mirgene){
         parmat <- parmat.list[[mirgene]]
       
         tparmat <- t(parmat)
         priormean <- pars$interactions$intercoeffs %*% tparmat
         ## priormeanVAR <- sapply(info$genes$m,function(mgene){
         ##   sum(diag( t(parmat[mgene,,drop=FALSE]) %*%
         ##            parmat[mgene,,drop=FALSE] %*%
         ##            intercoeffPrecInv ))
         ## })
         priormeanVAR <- sapply(info$genes$m,function(mgene){
           sum(diag( tparmat[,mgene,drop=FALSE] %*%
                    parmat[mgene,,drop=FALSE] %*%
                    intercoeffPrecInv ))
         })
         ret <- -0.5*log(2*pi) +
           0.5 * log(info$prediction.weight) +
             0.5 * gamExpLogx(pars$interactions$prec.common$alf,
                              pars$interactions$prec.common$bta) -
                                0.5 * info$prediction.weight *
                                  expect$interactions$prec.common *
                                    ( (pars$interactions$individual[,mirgene] -
                                       priormean)^2 +
                                     1/pars$interactions$individual.prec[,mirgene] +
                                     priormeanVAR )

         ## ret2 <- sapply(info$genes$m,function(mgene){
         ##   ## mirname <- info$possmir.genes$mir[mirgene]
         ##   ## mname <- info$possmir.gene$m[mgene]
         ##   ## parmat <- getParmat(mname,mirname,possmirs,info)
         ## 
         ##   parmat <- parmat.list[[mirgene]][mgene,,drop=FALSE]
         ##   
         ##   priormean <- sum( t(pars$interactions$intercoeffs) %*% t(parmat) )
         ##   priormeanVAR <- sum(diag( t(parmat) %*% parmat %*%
         ##                            intercoeffPrecInv ))
         ##   ret <- -0.5*log(2*pi) +
         ##     0.5 * log(info$prediction.weight) +
         ##       0.5 * gamExpLogx(pars$interactions$prec.common$alf,
         ##                        pars$interactions$prec.common$bta) -
         ##                          0.5 * info$prediction.weight *
         ##                            expect$interactions$prec.common *
         ##                              ( (pars$interactions$individual[mgene,mirgene] -
         ##                                 priormean)^2 +
         ##                               1/pars$interactions$individual.prec[mgene,mirgene] +
         ##                               priormeanVAR )
         ##   return(ret)
         ## })
         return(ret)
       }),
       
       qlogq=
       gaussEnt(pars$interactions$individual.prec)
       )
}

calcllh.interactions.intercoeffs <- function(priors,pars,info){
  expect <- getExpectations( pars )
  intercoeffPrecInv <- ginv(pars$interactions$intercoeffs.prec)

  list(
       qlogp=
       -0.5*length(info$interactionfactors)*log(2*pi) +
       0.5*sum(diag(as.matrix(log(priors$interactions$intercoeffs$prec)))) - 
       0.5 * ( t( pars$interactions$intercoeffs -
                 priors$interactions$intercoeffs$mean ) %*%
              priors$interactions$intercoeffs$prec %*% 
              as.matrix( pars$interactions$intercoeffs -
                        priors$interactions$intercoeffs$mean ) +
              sum(diag(as.matrix(priors$interactions$intercoeffs$prec %*%
                       intercoeffPrecInv) )) ),
       
       qlogq=
       multiGaussEnt(pars$interactions$intercoeffs.prec)
       )
}

calcllh.interactions.prec.common <- function(priors,pars,info){
  expect <- getExpectations( pars )
  list(
       qlogp=
       priors$interactions$prec.common$alf *
       log(priors$interactions$prec.common$bta) -
       lgamma( priors$interactions$prec.common$alf ) +
       ( priors$interactions$prec.common$alf - 1 ) *
       gamExpLogx(pars$interactions$prec.common$alf,
                  pars$interactions$prec.common$bta) -
       priors$interactions$prec.common$bta *
       expect$interactions$prec.common,
       
       qlogq=
       gamEnt(pars$interactions$prec.common$alf,
              pars$interactions$prec.common$bta)
       )
}


## calcllh.interactions.prec.individual <- function(priors,pars,info){
##   expect <- getExpectations( pars )
##   list(
##        qlogp=
##        priors$interactions$prec.individual$alf *
##        log(priors$interactions$prec.individual$bta) -
##        lgamma( priors$interactions$prec.individual$alf ) +
##        ( priors$interactions$prec.individual$alf - 1 ) *
##        gamExpLogx(pars$interactions$prec.individual$alf,
##                   pars$interactions$prec.individual$bta) -
##        priors$interactions$prec.individual$bta *
##        expect$interactions$prec.individual,
##        
##        qlogq=
##        gamEnt(pars$interactions$prec.individual$alf,
##               pars$interactions$prec.individual$bta)
##        )
## }
## 
## calcllh.interactions.prec.cluster <- function(priors,pars,info){
##   expect <- getExpectations( pars )
##   list(
##     qlogp=
##     priors$interactions$prec.cluster$alf *
##     log(priors$interactions$prec.cluster$bta) -
##     lgamma( priors$interactions$prec.cluster$alf ) +
##     ( priors$interactions$prec.cluster$alf - 1 ) *
##     gamExpLogx(pars$interactions$prec.cluster$alf,
##                pars$interactions$prec.cluster$bta) -
##     priors$interactions$prec.cluster$bta *
##     expect$interactions$prec.cluster,
## 
##     qlogq=
##     gamEnt(pars$interactions$prec.cluster$alf,
##             pars$interactions$prec.cluster$bta)
##     )
## }


calcllh.techprec <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  list(
       qlogp=
       lapply(list(spot="spot",slide="slide"),function(x){
         ( priors[[side]]$techprec[[x]]$alf * log(priors[[side]]$techprec[[x]]$bta) -
          lgamma( priors[[side]]$techprec[[x]]$alf ) +
          ( priors[[side]]$techprec[[x]]$alf - 1 ) *
          gamExpLogx(pars[[side]]$techprec[[x]]$alf,
                     pars[[side]]$techprec[[x]]$bta) -
          priors[[side]]$techprec[[x]]$bta *
          expect[[side]]$techprec[[x]] )
    }),
    
    qlogq=
    list(spot= gamEnt(pars[[side]]$techprec$spot$alf,
                       pars[[side]]$techprec$spot$bta),
         slide= gamEnt(pars[[side]]$techprec$slide$alf,
                        pars[[side]]$techprec$slide$bta) )
    )
}

  calcllh.clusprec <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  list(
    qlogp=
       priors[[side]]$clusprec$alf *
       gamExpLogx(pars[[side]]$clusprecBTAhyp$alf,
               pars[[side]]$clusprecBTAhyp$bta) -
    lgamma( priors[[side]]$clusprec$alf ) +
    ( priors[[side]]$clusprec$alf - 1 ) *
    gamExpLogx(pars[[side]]$clusprec$alf,
               pars[[side]]$clusprec$bta) -
    expect[[side]]$clusprecBTAhyp *
    expect[[side]]$clusprec,

    qlogq=
    gamEnt( pars[[side]]$clusprec$alf,
            pars[[side]]$clusprec$bta)
    )
}

calcllh.clusprecBTAhyp <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  list(
    qlogp=
    priors[[side]]$clusprecBTAhyp$alf *
    log(priors[[side]]$clusprecBTAhyp$bta) -
    lgamma( priors[[side]]$clusprecBTAhyp$alf ) +
    ( priors[[side]]$clusprecBTAhyp$alf - 1 ) *
    gamExpLogx(pars[[side]]$clusprecBTAhyp$alf,
               pars[[side]]$clusprecBTAhyp$bta) -
    priors[[side]]$clusprecBTAhyp$bta *
    expect[[side]]$clusprecBTAhyp,

    qlogq=
    gamEnt(pars[[side]]$clusprecBTAhyp$alf,
            pars[[side]]$clusprecBTAhyp$bta)
    )
}

calcllh.stageprec <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  if( any(unlist(info$typefollows[[side]]) %in% info$types) ){
    ret <- list(
      qlogp=
      priors[[side]]$stageprec$alf *
      gamExpLogx(pars[[side]]$stageprecBTAhyp$alf,
                 pars[[side]]$stageprecBTAhyp$bta) -
      lgamma( priors[[side]]$stageprec$alf ) +
      ( priors[[side]]$stageprec$alf - 1 ) *
      gamExpLogx(pars[[side]]$stageprec$alf,
                 pars[[side]]$stageprec$bta) -
      expect[[side]]$stageprecBTAhyp *
      expect[[side]]$stageprec,
      
      qlogq=
      gamEnt(pars[[side]]$stageprec$alf,
             pars[[side]]$stageprec$bta)
      )
  }else{
    ret <- list( qlogp=0, qlogq=0 )
  }
  return(ret)
}

calcllh.stageprecBTAhyp <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )

  if( any(unlist(info$typefollows[[side]]) %in% info$types) ){
    ret <- list(
      qlogp=
      priors[[side]]$stageprecBTAhyp$alf *
      log(priors[[side]]$stageprecBTAhyp$bta) -
      lgamma( priors[[side]]$stageprecBTAhyp$alf ) +
      ( priors[[side]]$stageprecBTAhyp$alf - 1 ) *
      gamExpLogx(pars[[side]]$stageprecBTAhyp$alf,
                 pars[[side]]$stageprecBTAhyp$bta) -
      priors[[side]]$stageprecBTAhyp$bta *
      expect[[side]]$stageprecBTAhyp,

      qlogq=
      gamEnt(pars[[side]]$stageprecBTAhyp$alf,
             pars[[side]]$stageprecBTAhyp$bta)
      )
  }else{
    ret <- list( qlogp=0, qlogq=0 )
  }
  return(ret)
}

calcllh.development <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  ## *development* isn't fully Bayesian (and isn't included
  ## in QlogQ below because it is optimized and not treated
  ## in the VB part of the algorithm. But, we leave it here
  ## in QlogP because then this function gives the (lower
  ## bound of) the log the marginal likelihood TIMES the
  ## likelihood of *development* according to its prior.           
  list(
    qlogp=
    unlist(
      lapply(names(pars[[side]]$development),function(typ){
        lapply(names(pars[[side]]$development[[typ]]),function(subtyp){
          ( priors[[side]]$development$alf *
           log(priors[[side]]$development$bta) -
           lgamma( priors[[side]]$development$alf ) +
           ( priors[[side]]$development$alf - 1 ) *
           log( abs(pars[[side]]$development[[typ]][subtyp]) ) -
           priors[[side]]$development$bta *
           abs(pars[[side]]$development[[typ]][[subtyp]]) )
        })
      }) ),

    qlogq= 0
    ## development= lapply(pars[[side]]$development,function(x){
    ##   sapply(x,function(y) -gamEnt(y$alf,y$bta) )
    ## }),
    )
}


calcllh.stagemean <- function(priors,pars,info,side=NULL,
                              intClusPrecInv=NULL){
  expect <- getExpectations( pars )
  list(
       qlogp=
       sapply( info$types,function(typ){
         if( "prior" %in% info$typefollows[[side]][[typ]] ){
           ret <- -0.5*log(2*pi) +
             0.5*log(priors[[side]]$stagemean$prec) - 
               0.5*priors[[side]]$stagemean$prec *
                 ( pars[[side]]$stagemean[,typ]^2 + 1/pars[[side]]$stagemean.prec[,typ] +
                  priors[[side]]$stagemean$mean^2 -
                  2 * pars[[side]]$stagemean[,typ]*priors[[side]]$stagemean$mean )
         }else{
           prevstagemean <- updateStageprec( xpr, pars, priors, info,
                                            side=side, forllh=TRUE,
                                            intClusPrecInv=intClusPrecInv)
           ret <- -0.5*log(2*pi) +
             0.5 * ( gamExpLogx(pars[[side]]$stageprec$alf,pars[[side]]$stageprec$bta) +
                    log(prevstagemean$allWTs[,typ]) ) -
                      0.5*expect[[side]]$stageprec * prevstagemean$allWTs[,typ] *
                        ( pars[[side]]$stagemean[,typ]^2 + 1/pars[[side]]$stagemean.prec[,typ] +
                         prevstagemean$psmeans[,typ]^2 + prevstagemean$psmean.vars[,typ] -
                         2 * pars[[side]]$stagemean[,typ]*prevstagemean$psmeans[,typ] )
         }
         return(ret)
       }),

       qlogq=
       gaussEnt(pars[[side]]$stagemean.prec)
       )
}

calcllh.slidemean <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  list(
    qlogp=
    lapply(info$types,function(typ){
      sapply(1:info$nclus[[side]],function(xi){
        { -0.5*log(2*pi) +
            0.5*gamExpLogx(pars[[side]]$clusprec$alf[xi,typ],
                           pars[[side]]$clusprec$bta[xi,typ]) -
                             0.5*expect[[side]]$clusprec[xi,typ] *
                               ( pars[[side]]$stagemean[xi,typ]^2 + 1/pars[[side]]$stagemean.prec[xi,typ] +
                                pars[[side]]$slidemean[,typ]^2 + 1/pars[[side]]$slidemean.prec[,typ] -
                                2 * pars[[side]]$stagemean[xi,typ]*pars[[side]]$slidemean[,typ] )
        }*pars[[side]]$clusmship[,xi]
      })
    }),

    qlogq=
    gaussEnt(pars[[side]]$slidemean.prec)
    )
}

calcllh.spotmean <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  list(
    qlogp=
    -0.5*log(2*pi) + 0.5*gamExpLogx(pars[[side]]$techprec$slide$alf,
                                    pars[[side]]$techprec$slide$bta) -
    0.5*expect[[side]]$techprec$slide *
    ( pars[[side]]$slidemean[,xpr[[side]]$targets$type]^2 +
     1/pars[[side]]$slidemean.prec[,xpr[[side]]$targets$type] +
     pars[[side]]$spotmean^2 + 1/pars[[side]]$spotmean.prec -
     2 * pars[[side]]$slidemean[,xpr[[side]]$targets$type]*pars[[side]]$spotmean ),

    qlogq=
    gaussEnt(pars[[side]]$spotmean.prec)
    )
}

calcllh.data <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  ## data=
  ## lapply(info$genes[[side]],function(gn){
  ##   -0.5*log(2*pi) +
  ##     0.5*gamExpLogx(pars[[side]]$techprec$spot$alf,
  ##                    pars[[side]]$techprec$spot$bta) -
  ##                      0.5*expect[[side]]$techprec$spot *
  ##                        apply( xpr[[side]]$E[info$genes[[side]]==gn,],1,function(xrow)
  ##                               ( xrow^2 +
  ##                                pars[[side]]$spotmean[gn,]^2 + 1/pars[[side]]$spotmean.prec[gn,] -
  ##                                2 * xrow*pars[[side]]$spotmean[gn,] ) )
  ## }),
  list(
    qlogp=
    -0.5*log(2*pi) + 0.5*gamExpLogx(pars[[side]]$techprec$spot$alf,
                                    pars[[side]]$techprec$spot$bta) -
    0.5*expect[[side]]$techprec$spot *
    ( xpr[[side]]$E^2 +
     pars[[side]]$spotmean[xpr[[side]]$genes[[info$gname[[side]]]],]^2 +
     1/pars[[side]]$spotmean.prec[xpr[[side]]$genes[[info$gname[[side]]]],] -
     2 * xpr[[side]]$E*pars[[side]]$spotmean[xpr[[side]]$genes[[info$gname[[side]]]],]),
    qlogq= 0
    )
}

calcllh.trend <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  list(
    qlogp=
    -0.5*log(2*pi) + 0.5*log(priors[[side]]$trend$prec) -
    0.5*priors[[side]]$trend$prec *
    ( pars[[side]]$trend^2 + 1/pars[[side]]$trend.prec +
     priors[[side]]$trend$mean^2 -
     2 * pars[[side]]$trend*priors[[side]]$trend$mean ),

    qlogq=
    ifelse( side=="mir",
           gaussEnt(pars[[side]]$trend.prec),0)
    )
}

calcllh.clusmship <- function(priors,pars,info,side=NULL){
  expect <- getExpectations( pars )
  list(
    qlogp=
    ifelse( priors[[side]]$clusmship==0,0,
           pars[[side]]$clusmship*log(priors[[side]]$clusmship) ),
    
    qlogq=
    -ifelse( pars[[side]]$clusmship==0,0,
           pars[[side]]$clusmship*log(pars[[side]]$clusmship) )
    )
}


calcllh <- function(priors,pars,info,possmirs){
  expect <- getExpectations( pars )
  llhlist <- NULL

  ## pre-calculation
  ##intercoeffPrecInv <- ginv(pars$interactions$intercoeffs.prec)
  
  llhlist$interactions <-
    list(
         cluster=calcllh.interactions.cluster(priors,pars,info),
         individual=calcllh.interactions.individual(priors,pars,info),
         intercoeffs=calcllh.interactions.intercoeffs(priors,pars,info),
         ## prec.individual=calcllh.interactions.prec.individual(priors,pars,info),
         ## prec.cluster=calcllh.interactions.prec.cluster(priors,pars,info)
         prec.common=calcllh.interactions.prec.common(priors,pars,info)
      )
    
  llhlist$mir <- NULL
  llhlist$m <- NULL
  for( side in c("mir","m") ){
    llhlist[[side]] <-
      list(
        techprec=calcllh.techprec(priors,pars,info,side=side),
        clusprec=calcllh.clusprec(priors,pars,info,side=side),
        clusprecBTAhyp=calcllh.clusprecBTAhyp(priors,pars,info,side=side),
        stageprec=calcllh.stageprec(priors,pars,info,side=side),
        stageprecBTAhyp=calcllh.stageprecBTAhyp(priors,pars,info,side=side),
        development=calcllh.development(priors,pars,info,side=side),
        stagemean=calcllh.stagemean(priors,pars,info,side=side),
        slidemean=calcllh.slidemean(priors,pars,info,side=side),
        spotmean=calcllh.spotmean(priors,pars,info,side=side),
        data=calcllh.data(priors,pars,info,side=side),
        ##trend=calcllh.trend(priors,pars,info,side=side),
        clusmship=calcllh.clusmship(priors,pars,info,side=side)
        )
  }

 
  ##
  ##llh <- sum(unlist(qlogp)) - sum(unlist(qlogqentropy))
  llh <- sum(unlist(llhlist))


  ## llh01 <- calcllhOLD(priors,pars,info,possmirs)
  ## cat(llh[[1]]-llh01[[1]],"\n")

  ## return(llh)
  return( list(llh,llhlist) )
}  

