## This function is intended to run many simulations and subsequent
## model fits.

clusteringRunFunction <- function(xpr=NULL,
                                  info=NULL,
                                  doclustering=FALSE,
                                  possmirs=NULL,
                                  reps=1){

 
  initialiters <- 5
  dataiters <- 1
  smalldataiters <- 3
  interactioniters <- 10
  bigiters <- 50

  ##doclustering <- TRUE

  ## initialiters <- 1
  ## dataiters <- 2
  ## smalldataiters <- 2
  ## interactioniters <- 5
  ## bigiters <- 1
  ## doclustering <- TRUE


  ## ###########################
  ## for many cores, make these equal, or multiples
  mc.cores <- 3
  simulruns <- 8

  stopruns <- 1


  clusrange.total <- list()
  clusrange.total$mir <- c(11,length(info$genes$mir))
  clusrange.total$m <- c(2,length(info$genes$m))


  for( side in c("mir","m")){
    if( !doclustering[[side]] ){
      clusrange.total[[side]] <- rep(length(info$genes[[side]]),2)
    }
  }


  ## #####################

  clusrange <- list()
  clusrange$mir <- clusrange.total$mir
  clusrange$m <- clusrange.total$m


  numruns.mat <- matrix(0,clusrange.total$mir[2],
                        clusrange.total$m[2])



  numruns.atoptim <- 0


  bestresult <- list(ret=NULL,
                     llh=-Inf,
                     nclus=lapply(info$genes,length))
    
  doside <- "mir"

  ## ## ######################################
  ## ## this is the big run; it can take a long time
  ## while( numruns.atoptim < stopruns ){

  otherside <- ifelse(doside=="mir","m","mir")

  doclusnums <- lapply( c(mir="mir",m="m"),function(side){
    if( side==doside ){
      y <- round( seq(from=clusrange[[side]][1],
                      to=clusrange[[side]][2],
                      length.out=simulruns) )
    }else{
      y <- bestresult$nclus[[side]]
    }
    return(y)
  })

  timestamp()
  cat( "Starting runs with nclus: mir--", doclusnums$mir,
      "; m--", doclusnums$m,"\n")

  ##clear the diagnostic file for the next run
  sink("fitModel_diagnostic.txt",append=FALSE)
  cat("restart\n")
  sink(); sink()

  info$nclus <- bestresult$nclus
  
  allret <- mclapply(doclusnums[[doside]],function(nc){
    ##nc <- 3
    info$nclus[[doside]] <- nc
    sink("fitModel_diagnostic.txt",append=TRUE)
    
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
      
    sink()

    return(ret2)
  },mc.cores=mc.cores)
  
  finalret <- list( pars=list(xpr=xpr,
                      info=info,
                      doclustering=doclustering,
                      possmirs=possmirs,
                      reps=reps),
                   allret=allret )
  return(finalret)
}
