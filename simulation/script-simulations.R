library(multicore)

## #####################################
## #####################################
## do a big lapply() for simulations

hugeret <- lapply( 1:5,function(x){

  
  ## #####################################
  ## file of the main user-defined parameters
  source("README_userDefinedParameters.R")

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


  source("clusteringRunFunction.R")

  reps <- 5
  
  runret <- list()

  info$interactionfactors <- c(constant="constant")
  info$prediction.weight <- 1
  runret$halfclus1 <- clusteringRunFunction(xpr=xpr,info=info,
                                           doclustering=list(mir=TRUE,m=FALSE),
                                           possmirs=NULL,reps=reps)

  info$interactionfactors <- c(constant="constant",pred="pred")
  info$prediction.weight <- 1
  runret$halfclus.pred1 <- clusteringRunFunction(xpr=xpr,info=info,
                                                doclustering=list(mir=TRUE,m=FALSE),
                                                possmirs=possmirs,reps=reps)

  hret <- list(simret=simret,
               runret=runret)
  return(hret)
})



## save(hugeret,file="five10x10nets_predwt4.Rdata")

hugeret.orig <- hugeret


## lapply(hugeret,function(x){
##   sapply(x$runret,function(y)
##          unlist(y$aucs[c("auroc","aupr")]))
## })

load("five10x10nets_predwt4.Rdata")


## ################################
## check accuracy
     

## aucs: all runs
all.ret <- lapply( hugeret, function(hr){
  lapply( hr$runret,function(rr){
    aucs.list <- lapply( rr$allret,function(xclusrun){
      sapply(xclusrun,function(yrep){
        xpars <- yrep$pars
        inferredints <- xpars$interactions$individual *
          xpars$interactions$individual.prec^0.5
        goldstd <- hr$simret$intermat
        aucs <- checkInteractionResults( inferredints, goldstd )
        return(aucs[c("auroc","aupr")])
      })
    })
  })
})


## aucs: aggregated, one per clusnum
best.ret <- lapply( hugeret, function(hr){
  lapply( hr$runret,function(rr){
    aucs.list <- sapply( rr$allret,function(xclusrun){

      bestllh <- which.max(sapply(xclusrun,function(yrep) yrep$llh[[1]] ))
      xpars <- xclusrun[[bestllh]]$pars
      inferredints <- xpars$interactions$individual *
        xpars$interactions$individual.prec^0.5

      ## ## average of rankings (probably not good)
      ## intsum <- matrix(0,nrow(x[[1]]$pars$interactions$individual),
      ##                  ncol(x[[1]]$pars$interactions$individual))
      ## intvarsum <- intsum
      ## for( i in 1:length(x) ){
      ##   intsum <- intsum + x[[i]]$pars$interactions$individual
      ##   intvarsum <- intvarsum + x[[i]]$pars$interactions$individual.prec^-1
      ## }
      ## inferredints <- intsum / intvarsum^0.5
      
      goldstd <- hr$simret$intermat
      aucs <- checkInteractionResults( inferredints, goldstd )
      return(aucs[c("auroc","aupr")])
    })
  })
})


best.ret.auroc <- lapply(best.ret,function(br){
  t( sapply(br,function(x) round(unlist(x[1,]),digits=3)) )
})

## names(best.ret[[1]])
use.fit.ind <- 2
br.auroc.simple <- t( sapply(best.ret,function(br){
  round(unlist(br[[use.fit.ind]][1,]),digits=3) 
}) )


for( i in 1:length(best.ret.auroc) ){
  x11()
  matplot(x=1:8,y=t(best.ret.auroc[[i]]),type="l")
}



## ####################################
## proportions over/under-performed, all
overunder <- lapply(all.ret,function(xnet){
  lapply(xnet,function(xpar){
    sapply(xpar,function(xcn){
      c(auroc.gt=mean( unlist(xcn["auroc",])>unlist(xpar[[length(xpar)]]["auroc",])),
        auroc.lt=mean( unlist(xcn["auroc",])<unlist(xpar[[length(xpar)]]["auroc",])) ##,
        ##aupr=mean(unlist(xcn["aupr",]>=unlist(xpar[[length(xpar)]]["aupr",])))
        )
    })
  })
})


## same as above, rearranged
runnames <- names(all.ret[[1]])
names(runnames) <- runnames
perf.tabs.gt <- lapply(runnames,function(rn){
  t( sapply(all.ret,function(xnet){
    xpar <- xnet[[rn]]
    sapply(xpar,function(xcn){
      c(auroc.gt=mean( unlist(xcn["auroc",])>unlist(xpar[[length(xpar)]]["auroc",])) ##,
        ##auroc.lt=mean( unlist(xcn["auroc",])<unlist(xpar[[length(xpar)]]["auroc",])) ##,
        ##aupr=mean(unlist(xcn["aupr",]>=unlist(xpar[[length(xpar)]]["aupr",])))
        )
    })
  }) )
})
perf.tabs.lt <- lapply(runnames,function(rn){
  t( sapply(all.ret,function(xnet){
    xpar <- xnet[[rn]]
    sapply(xpar,function(xcn){
      c(##auroc.gt=mean( unlist(xcn["auroc",])>unlist(xpar[[length(xpar)]]["auroc",])) ##,
        auroc.lt=mean( unlist(xcn["auroc",])<unlist(xpar[[length(xpar)]]["auroc",])) ##,
        ##aupr=mean(unlist(xcn["aupr",]>=unlist(xpar[[length(xpar)]]["aupr",])))
        )
    })
  }) )
})

t(sapply(perf.tabs.gt,function(x)
       apply(x,2,mean)))[use.fit.ind,]
t(sapply(perf.tabs.lt,function(x)
       apply(x,2,mean)))[use.fit.ind,]
       


## #######################################
## tallying performance of the best among reps

lapply(best.ret,function(xnet){
  lapply(xnet,function(xpar){
    rbind(auroc=ifelse(unlist(xpar["auroc",])>unlist(xpar["auroc",ncol(xpar)]),1,0) ##,
      ##aupr=ifelse(unlist(xpar["aupr",]>=unlist(xpar["aupr",ncol(xpar)])),1,0))
          )
  })
})


## same as above, rearranged
runnames <- names(best.ret[[1]])
names(runnames) <- runnames
perf.tabs.gt <- lapply(runnames,function(rn){
  t( sapply(best.ret,function(xnet){
    xpar <- xnet[[rn]]
    c(auroc.gt=ifelse(unlist(xpar["auroc",])>unlist(xpar["auroc",ncol(xpar)]),1,0) )
  }) )
})
perf.tabs.lt <- lapply(runnames,function(rn){
  t( sapply(best.ret,function(xnet){
    xpar <- xnet[[rn]]
    c(auroc.lt=ifelse(unlist(xpar["auroc",])<unlist(xpar["auroc",ncol(xpar)]),1,0) )
  }) )
})

t(sapply(perf.tabs.gt,function(x)
       apply(x,2,mean)))[use.fit.ind,]
t(sapply(perf.tabs.lt,function(x)
       apply(x,2,mean)))[use.fit.ind,]



## ########################################
## see where predictions helped


runnames <- names(best.ret[[1]])
names(runnames) <- runnames

perf.tabs.pred <- lapply(best.ret,function(xnet){
  t( sapply(1:(length(runnames)/2),function(rn){
    xpar <- xnet[[2*rn-1]]
    xpar.pred <- xnet[[2*rn]]
    c(auroc.gt=ifelse(unlist(xpar.pred["auroc",])>unlist(xpar["auroc",ncol(xpar)]),1,0) )
  }) )
})

sapply(perf.tabs.pred,mean)





## apply(xpr$mir$E,1,function(x) acf(x,lag.max=1,plot=FALSE))


## ################################################
## make a table of the probability of targeting and then sort it

source(paste(sourcedir,"generateNiceResults.R",sep="/"))

tablelist <- generateNiceResults(xpr,llh,
                                 pars,priors,
                                 info,
                                 possmirs)

write.csv(tablelist$intz,file="resultsTable.csv")




## ## ################################################
## ## ################################################
## ## ################################################
## ## ################################################
## ## mirname <- "mmu-miR-200b"
## ## mname <- "Zeb2"
## ## mirgene <- names(info$possmir.genes$mir[info$possmir.genes$mir==mirname])
## ## mgene <- names(info$possmir.genes$m[info$possmir.genes$m==mname])
## ## checkclusts <- list()
## ## checkclusts$mir <- which(pars$mir$clusmship[mirgene,]>0.9)
## ## checkclusts$m <- which(pars$m$clusmship[mgene,]>0.9)


for( i in 1:3 ){

  ## ret <- allret[[i]]
  ret <- hugeret[[1]]$runret$halfclus[[i]][[1]]


  
  pars <- ret$pars
  priors <- ret$priors
  info <- ret$info
  llh <- ret$llh

  
  expect <- getExpectations( pars )
  ## expect[[side]]$clusprec

  checkclusts <- list()
  checkclusts$mir <- 1:info$nclus$mir
  checkclusts$m <- 1:info$nclus$m

  x11()
  par(mfrow=c(2,1),mar=c(2,2,1,1))
  for( side in c("mir","m") ){
    ##doclusts <- 1:info$nclus[[side]]
    ##doclusts <- 2
    doclusts <- checkclusts[[side]]
    matplot( t(pars[[side]]$stagemean[doclusts,,drop=FALSE]),type="l",lwd=4,lty=2,
            col=doclusts)
    nulldev <- lapply(doclusts,function(i){
      clusmembs <- info$genes[[side]][pars[[side]]$clusmship[,i]>0.5]
      if( length(clusmembs)>0 ){
        geneprofs <- t(pars[[side]]$slidemean[clusmembs,,
                                              drop=FALSE])
        matlines(geneprofs,col=i,lty=2)
        sp2 <- sapply(info$types,function(typ){
          spotprofs <- pars[[side]]$spotmean[clusmembs,
                                             xpr[[side]]$targets$type==typ,
                                             drop=FALSE]
          xval <- which(info$types==typ) +
            0.1 * ( col(spotprofs) - mean(col(spotprofs)) )
          points(xval,spotprofs,col=i,lty=2,pch="o")
        })
        d2 <- sapply(info$types,function(typ){
          usespots <- which(xpr[[side]]$genes[[info$gname[[side]]]] %in%
                            clusmembs)
          datpts <- xpr[[side]]$E[usespots,
                                  xpr[[side]]$targets$type==typ,
                                  drop=FALSE]
          xval <- which(info$types==typ) +
            0.1 * ( col(datpts) - mean(col(datpts)) )
          points(xval,datpts,col=i,lty=2,pch="x")
        })
      }
    })
    matlines( t(pars[[side]]$stagemean[doclusts,,drop=FALSE]),type="l",lwd=4,lty=2,
             col=doclusts)
  }

  
}
