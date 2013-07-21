## ###################################################
## ###################################################
##
## This file contains user-specified parameters that
## the R script needs to run. Read the parameter
## descriptions carefully.
## 
## ###################################################
## ###################################################


## #######
## specify the directories where important things are
##
## this is where the src/ files are (modelFunctions.R, etc)
sourcedir <- "./src" # default location from git repo

## location of the files from TargetScan, miRanda, etc.
databasefiledir <- "./src/databasefiles" # also in git repo

## a place to look for and/or save files we generate (e.g. f-test
## results).  you should make this directory if it doesn't
## exist. Also, delete its contents if you change data sets. (Not
## necessary for multiple runs on the same starting data set, but with
## different parameters).
generatedfiledir <- "./generatedfiles"


## ####################
## choose organism and array platform, or set both to "simulation"
## for a single file with all expression values, set
## arrayplatform <- "expressiontable"

data.organism <- "simulation"
arrayplatform <- "simulation"

simnclus <- list(mir=25,
                 m=25)


## OR

## ######
## ## uncomment one organism and comment out the other
##data.organism <- "mouse"
##data.organism <- "human"

## ######
## ##uncomment one array platform and comment out the other
##arrayplatform <- "agilent"
##arrayplatform <- "affymetrix" # right now, affy works only with human
##arrayplatform <- "expressiontable"



## #############################

## specify data either from platform-specific files (no guarantee that
## they can be read) or from a single file (for each of miR and mRNA
## data) with a table of expression values complete one section or the
## other


## ##############
## PLATFORM SPECIFIC FILES:
## Tab-separated file with three columns with these headers:
## filename -- complete path to a data file
## RNAtype -- either "m" or "miR"
## type -- samples of the same type are considered replicates
## TableOfSamplesFilename <- "TableOfSamples.txt"
##
## column header in the data files that we want to use as main list of
## probes or genes; probe names tend to be better
##
## gname <- list(mir="SystematicName",m="GeneName")
## this has worked for Agilent array data
gname <- list(mir="ProbeName",
              m="ProbeName")
##
##
## column name in the data files containing gene/miR names that match
## those in the database (targetScan, miRanda, etc) files
##
## this has worked for Agilent array data
possmir.gname <- list(mir="SystematicName",
                      m="GeneName")



## ########
## SINGLE FILE WITH EXPRESSION TABLE FOR EACH OF mRNA and miR:
## make sure arrayplatform <- "expressiontable" above
##
## tab-separated file with samples in columns, probes in rows
## column names (first row) must match sample/stage types below
## row names should not be included
exprmat.filename.mir <- "exprmatFile_mir.txt"
exprmat.filename.m <- "exprmatFile_m.txt"
##
## This file gives the probe names (first column) and gene
## names(second column), in the same order they appear in the
## expression table. Gene names should match those in the databases
## (e.g. TargetScan)
probenames.filename.mir <- "probenamesFile_mir.txt"
probenames.filename.m <- "probenamesFile_m.txt"


if( arrayplatform=="expressiontable" ){
  gname <- list(mir="ProbeName",
                m="ProbeName")
  possmir.gname <- list(mir="GeneName",
                        m="GeneName")
}

## use this line if you want to make the expression table files for
## your data AFTER you've already loaded the data once and *xpr.orig*
## exists
##
##source(paste(sourcedir,"makeExpressionTableFile.R",sep="/"))



## ###############################
## ordered list of sample types
##

## Every sample "type" (from the user-defined table of samples) must
## either precede or follow some other type. Make sure there are no
## cycles (e.g. don't let something precede something else which
## precedes the original type). Also, every type must follow
## something; use "prior" for initial types. Which stages "precede"
## and "follow" aren't nearly as important as the connections between
## types. The goal is to make connections between types that we know
## should be similar. Multiple parents and children are allowed. Do
## not use special characters, spaces, or ".", or anything but lettters and
## numbers in the types. MUST START WITH A LETTER.

## use the format
## [type] = "type it follows"
##
## for example:


typefollows <- list(s01="prior",
                    s02="s01",
                    s03="s02",
                    s04="s03",
                    s05="s04",
                    s06="s05",
                    s07="s06",
                    s08="s07" #,
                    ## s09="s08",
                    ## s10="s09"
                    )

## typefollows <- list(healthyDay1="prior",
##                     healthyDay2="healthyDay1",
##                     inductionType1Day3="healthyDay2",
##                     inductionType1Day4="inductionType1Day3",
##                     inductionType2Day3="healthyDay2",
##                     inductionType2Day4=c("inductionType2Day3","inductionType2Day3"))

## typefollows <- list(LTHSC="prior",
##                     STHSC="LTHSC",
##                     MPP=c("STHSC","KSL"),
##                     CMP="MPP",
##                     MEP="CMP",
##                     GMP="CMP",
##                     ## CLP="MPP",
##                     KSL=c("LTHSC","STHSC") )

##   OR    ##


## instead of specifying *typefollows* above, you can choose a
## reference type, and then *typefollows* will be generated with all
## types following that reference type

## typefollows <- NULL
## reference <- "LTHSC"


## ####################
##
## Threshold for limiting the data set size to the genes/miRs whose
## expression changes significantly between types. Running the code
## takes a very long time if the number of miRs/genes is much larger
## than ~200, but most expression profiles don't change significantly
## between types, so it's OK to throw out the profiles that don't
## change. They won't affect results anyway. A good starting point is
## to use ftestthresh=0.05 and padjustmethod="none", but that is very
## generous. If there are still lots of miRs (>>100) or genes (>>300),
## switch to "fdr" and then gradually decrease *ftestthresh* until you
## get reasonable numbers. If you have a lot of types, *ftestthresh*
## can be very small, e.g. 1e-6 for 8 types with 2-3 replicates each.

##ftestthresh <- 1e-5

ftestthresh <- list()
ftestthresh$mir <- 100
ftestthresh$m <- 300

##padjustmethod <- "fdr"
padjustmethod <- "none"



## ###############
## include only miRs and mRNAs predicted to be involved in targeting?
##
## choosing TRUE dramatically reduces the data set, but puts more
## emphasis on the databases

predictedonly <- FALSE

prediction.weight <- 1




