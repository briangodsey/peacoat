individual <- TRUE
sequential <- TRUE


typefollows <- list( ##normal="prior",
                    IA="prior",
                    IIA="IA",
                    IIIA="IIA",
                    IIIB="IIIA",
                    PCL="IA")

## ############################################
## ############################################
## This changes the model design based on
## *individual* and *sequential* from above
## ############################################
## ############################################
if( individual ){

  types.indiv <- list()
  for( side in c(mir="mir",m="m") ){
    patnums <- sapply(as.character(xpr.orig[[side]]$targets$source_name_ch1),function(x){
      y <- strsplit(x," ")[[1]]
      return( y[length(y)] )
    })
    types.indiv[[side]] <- paste(types[[side]],patnums,sep="_")
  }

  matchedtypes <- unique(intersect(types.indiv$mir,
                                   types.indiv$m))
  names(matchedtypes) <- matchedtypes

  matchedtypes.sub <- sapply(matchedtypes,function(mt)
                             strsplit(mt,"_")[[1]][1] )

  ## get rid of types/individuals not in *typefollows* above
  ind <- which( matchedtypes.sub %in% names(typefollows) )
  matchedtypes <- matchedtypes[ind]
  matchedtypes.sub <- matchedtypes.sub[ind]

  if( sequential ){
    ## SEQUENTIAL DESIGN BY DISEASE STAGE AND PATIENT NUMBER
    typefollows.indiv <- lapply(matchedtypes,function(mt){
      bigtype <- strsplit(mt,"_")[[1]][1]
      tfoll <- typefollows[[bigtype]]
      if( tfoll=="prior" ){
        ret <- "prior"
      }else{
        ret <- matchedtypes[matchedtypes.sub==tfoll]
      }
      return(ret)
    })
  }else{ ## if NOT sequential
    ## IGNORE SEQUENCE BUT SEPARATE PATIENTS
    typefollows.indiv <- lapply(matchedtypes,function(mt){
      bigtype <- strsplit(mt,"_")[[1]][1]
      tfoll <- typefollows[[bigtype]]
      if( tfoll=="prior" ){
        ret <- "prior"
      }else{
        ret <- matchedtypes[matchedtypes.sub=="IA"]
      }
      return(ret)
    })
  }
    
  ## replace variables with the individual versions
  types <- types.indiv
  typefollows <- typefollows.indiv
  xpr.orig$mir$targets$type <- types.indiv$mir
  xpr.orig$m$targets$type <- types.indiv$m

}else{  ## if NOT individual
  if( sequential ){
    ## there are no "normal" mRNA arrays (DOH)
    typefollows <- list( ##normal="prior",
                        IA="prior",
                        IIA="IA",
                        IIIA="IIA",
                        IIIB="IIIA",
                        PCL="IA")
  }else{
    typefollows <- list( ##normal="prior",
                        IA="prior",
                        IIA="IA",
                        IIIA="IA",
                        IIIB="IA",
                        PCL="IA")
  }
}
