library(limma)


sampletable <- read.table(TableOfSamplesFilename,
                          sep="\t",header=TRUE,stringsAsFactor=FALSE)


fullfns.m <- sampletable$filename[sampletable$RNAtype=="m"]
fullfns.mir <- sampletable$filename[sampletable$RNAtype=="mir"]


## ####################################################
## loading of data

##allexpression <- read.maimages(files=fullfns[1],source="agilent")

## mRNA

headers.m <- c("FEATURES", "FeatureNum", "Row", "Col", "accessions",
               "chr_coord", "SubTypeMask", "SubTypeName", "Start",
               "Sequence", "ProbeUID", "ControlType", "ProbeName",
               "GeneName", "SystematicName", "Description",
               "PositionX", "PositionY", "gSurrogateUsed", "gIsFound",
               "gProcessedSignal", "gProcessedSigError",
               "gNumPixOLHi", "gNumPixOLLo", "gNumPix", "gMeanSignal",
               "gMedianSignal", "gPixSDev", "gPixNormIQR",
               "gBGNumPix", "gBGMeanSignal", "gBGMedianSignal",
               "gBGPixSDev", "gBGPixNormIQR", "gNumSatPix",
               "gIsSaturated", "gIsLowPMTScaledUp",
               "gIsFeatNonUnifOL", "gIsBGNonUnifOL", "gIsFeatPopnOL",
               "gIsBGPopnOL", "IsManualFlag", "gBGSubSignal",
               "gBGSubSigError", "gIsPosAndSignif", "gPValFeatEqBG",
               "gNumBGUsed", "gIsWellAboveBG", "gBGUsed", "gBGSDUsed",
               "ErrorModel", "gSpatialDetrendIsInFilteredSet",
               "gSpatialDetrendSurfaceValue", "SpotExtentX",
               "SpotExtentY", "gNetSignal", "gMultDetrendSignal",
               "gProcessedBackground", "gProcessedBkngError",
               "IsUsedBGAdjust", "gInterpolatedNegCtrlSub",
               "gIsInNegCtrlRange", "gIsUsedInMD")

## miRNA

headers.mir <- c("FEATURES", "FeatureNum", "Row", "Col",
                 "SubTypeMask", "ControlType", "ProbeName",
                 "SystematicName", "PositionX", "PositionY",
                 "gProcessedSignal", "gProcessedSigError",
                 "gNumPixOLHi", "gNumPixOLLo", "gNumPix",
                 "gMeanSignal", "gPixSDev", "gMedianSignal",
                 "gBGMedianSignal", "gBGPixSDev", "gIsSaturated",
                 "gIsLowPMTScaledUp", "gIsFeatNonUnifOL",
                 "gIsBGNonUnifOL", "gIsFeatPopnOL", "gIsBGPopnOL",
                 "IsManualFlag", "gBGSubSignal", "gIsPosAndSignif",
                 "gIsWellAboveBG", "SpotExtentX", "gBGMeanSignal",
                 "gTotalProbeSignal", "gTotalProbeError",
                 "gTotalGeneSignal", "gTotalGeneError",
                 "gIsGeneDetected")

ann <- unique(c('FeatureNum ','Row','Col','accessions',
                'chr_coord','SubTypeMask','SubTypeName',
                'ProbeUID','ControlType', 'ProbeName',
                'GeneName','SystematicName','accessions','chr_coord'))
ann.m <- intersect(ann,headers.m)
ann.mir <- intersect(ann,headers.mir)


## ##########################################
## load the data from the files -- TAKES TIME
xpress.m <- read.maimages(fullfns.m, source="generic",
                      columns=list(Rf="gMedianSignal", Rb="gBGMedianSignal",
                        ##Gf="gMeanSignal", Gb="gBGMeanSignal",
                        Flag="IsManualFlag",
                        Xcoord="PositionX",Ycoord="PositionY"),
                      annotation=ann.m)
xpress.m$targets <- sampletable[sampletable$RNAtype=="m",]


fnames <- sapply(xpress.m$targets$filename,function(fn){
  fnsplit <- strsplit(fn,split="/")[[1]]
  ret <- fnsplit[length(fnsplit)]
})
colnames(xpress.m$E) <- fnames
colnames(xpress.m$Eb) <- fnames


  

xpress.mir <- read.maimages(fullfns.mir, source="generic",
                        columns=list(Rf="gMedianSignal", Rb="gBGMedianSignal",
                          ##Gf="gMeanSignal", Gb="gBGMeanSignal",
                          Flag="IsManualFlag",
                          Xcoord="PositionX",Ycoord="PositionY"),
                        annotation=ann.mir)
xpress.mir$targets <- sampletable[sampletable$RNAtype=="mir",]

fnames <- sapply(xpress.mir$targets$filename,function(fn){
  fnsplit <- strsplit(fn,split="/")[[1]]
  ret <- fnsplit[length(fnsplit)]
})
colnames(xpress.mir$E) <- fnames
colnames(xpress.mir$Eb) <- fnames


## #################################################
## normalization

cat("\nQuantile normalizing...\n\n Don't worry about 'non-standard columns' warning. \n\n")

##xpress.m.scale <- normalizeBetweenArrays(xpress.m,method="scale")
xpress.m.qu <- normalizeBetweenArrays(xpress.m,method="quantile")

## xpress.mir.scale <- normalizeBetweenArrays(xpress.mir,method="scale")
xpress.mir.qu <- normalizeBetweenArrays(xpress.mir,method="quantile")



## choose expression objects to use as the main data, one for *mir*
## and one for *m*
xpr.orig <- NULL
xpr.orig$mir <- xpress.mir.qu
## xpr.orig$mir <- xpress.mir.scale
xpr.orig$m <- xpress.m.qu
## xpr.orig$m <- xpress.m.scale
