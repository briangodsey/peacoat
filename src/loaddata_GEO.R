## source("http://www.bioconductor.org/biocLite.R")
## biocLite("GEOquery")

library(limma)
library(Biobase)
library(GEOquery)


gse.obj <- lapply(geo.files,function(fn)
                  getGEO(filename=fn) )
platforms <- sapply(gse.obj,annotation)
gpl.obj <- lapply(platforms, function(plat)
                  getGEO(plat,destdir=".")) ##datafiledir))


xpress <- lapply(c(mir="mir",m="m"), function(side){
  ret <- list()
  ret$E <- log(exprs(gse.obj[[side]]))
  genestab <- as.data.frame( apply(Table(gpl.obj[[side]]),
                                   2,as.character),
                            stringsAsFactors=FALSE)
  rownames(genestab) <- genestab$ID
  ##indorder <- sapply(rownames(ret$E),function(x) which(genestab$ID==x))
  ret$genes <- genestab[rownames(ret$E),] ## genestab[indorder,]
  ret$targets <- pData(gse.obj[[side]])
  ret1 <- new("EList",ret)
  return(ret1)
})


## #################################################
## normalization

## quantile
xpress.qu <- lapply( xpress, function(x){
  ret <- x
  ret$E <- normalizeBetweenArrays(x$E,method="quantile")
  return(ret)
})

xpr.orig <- xpress.qu


