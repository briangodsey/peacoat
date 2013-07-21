types <- NULL


## ## this would copy a certain column into the "type" column
## types$mir <- as.character(xpr.orig$mir$targets[,"characteristics_ch1.3"])
## types$m <- as.character(xpr.orig$m$targets[,"characteristics_ch1.2"])


## here, we use the "durie-salmon" stage to order the data
## we leave out "nd"
types$mir <- sapply(as.character(xpr.orig$mir$targets[,"characteristics_ch1.3"]),
                    function(x){
                      if( x == "" ){
                        ret <- "normal"
                      }else{
                        ret0 <- gsub("stage (durie-salmon): ","",x,fixed=TRUE)
                        ret1 <- gsub("I A","IA",ret0,fixed=TRUE)
                        ret <- gsub("I B","IB",ret1,fixed=TRUE)
                      }
                      return(ret)
                    })
types$m <- sapply(as.character(xpr.orig$m$targets[,"characteristics_ch1.2"]),
                  function(x){
                    if( x == "" ){
                      ret <- "normal"
                    }else{
                      ret0 <- gsub("stage (durie-salmon): ","",x,fixed=TRUE)
                      ret1 <- gsub("I A","IA",ret0,fixed=TRUE)
                      ret <- gsub("I B","IB",ret1,fixed=TRUE)
                    }
                    return(ret)
                  })
names(types$mir) <- types$mir
names(types$m) <- types$m



## ###################
xpr.orig$mir$targets$type <- types$mir
xpr.orig$m$targets$type <- types$m

