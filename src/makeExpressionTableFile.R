
if( !exists("xpr.orig") ){
  stop("Data has not been loaded; xpr.orig does not exist.\n")
}

exprtab.list <- lapply( c(mir="mir",m="m"),function(side){
  xtab <- xpr.orig[[side]]$E
  rownames(xtab) <- xpr.orig[[side]]$genes[,info$gname[[side]]]
  colnames(xtab) <- xpr.orig[[side]]$targets$type
  return(xtab)
})

probenames.list <- lapply( c(mir="mir",m="m"),function(side){
  probenam <- xpr.orig[[side]]$genes[,c(info$gname[[side]],info$possmir.gname[[side]])]
  ##probenam <- probenam[!duplicated(probenam),]
  return(probenam)
})


for( side in c(mir="mir",m="m") ){
  write.table( exprtab.list[[side]],
              file=paste("exprmatFile_",side,".txt",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE)
  write.table( probenames.list[[side]],
              file=paste("probenamesFile_",side,".txt",sep=""),
              quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}
