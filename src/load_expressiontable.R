library(limma)

xprmat.fns <- NULL
probenam.fns <- NULL
xprmat.fns$mir <- exprmat.filename.mir
xprmat.fns$m <- exprmat.filename.m
probenam.fns$mir <- probenames.filename.mir
probenam.fns$m <- probenames.filename.m

expertab.list <- NULL
probenames.list <- NULL
for( side in c(mir="mir",m="m") ){
  expertab.list[[side]] <- read.table( xprmat.fns[[side]],
                                      as.is=TRUE,header=TRUE,sep="\t")
  probenames.list[[side]] <- read.table( probenam.fns[[side]],
                                      as.is=TRUE,header=FALSE,sep="\t",
                                      row.names=NULL)
}



cat("\nQuantile normalizing the data and building the main data object.\n")
xpr.orig <- list()
gname <- list()
possmir.gname <- list()
for( side in c(mir="mir",m="m") ){
  xpr.orig[[side]] <- NULL
  xpr.orig[[side]]$E <- normalizeBetweenArrays(as.matrix(expertab.list[[side]]),
                                               method="quantile")
  xpr.orig[[side]]$genes <- probenames.list[[side]]
  colnames(xpr.orig[[side]]$genes) <- c("ProbeName","GeneName")
  
  xpr.orig[[side]]$targets <- data.frame(sapply(colnames(expertab.list[[side]]),function(x){
    strsplit(x,split=".",fixed=TRUE)[[1]][1]
  }),stringsAsFactors=FALSE )
  colnames(xpr.orig[[side]]$targets) <- "type"

  xpr.orig[[side]] <- new("EList",xpr.orig[[side]])

  gname[[side]] <- colnames(xpr.orig[[side]]$genes)[1]
  possmir.gname[[side]] <- colnames(xpr.orig[[side]]$genes)[2]

}



