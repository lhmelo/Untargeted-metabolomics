#title: "XCMS.function"
#author: "Laura Hmelo"
#last update: "January 27, 2017"
#
#This is a version of XCMS.function.R which merges in features of Katherine Heal's XCMS function.  Namely, it brings in capability to run XCMS on all four metabolite fractions in a looped fraction.  It removes parameter definitions from heading of Example_workflow to Params.csv.

# 
# xcmsfunction----
# This script uses XCMS to generate a list of peaks that have been
# aligned, RT corrected, realigned, and recursively filled. 
# 
# Before running the function:
# 1. Define ResultsDIR and DataDIR
# 2. Define extractiontype in function call
# 
# Input:
# 1. .mzXML sample files unbinned in a single folder
# 
# Output is a list containing:
# 1. a processed xcmsSet object (xset3)
# 2. a data.frame containing a list of peaks that have been
# aligned, RT corrected, realigned, and recursively filled (xset.allpeaks)

xcmsfunction <- function() {

library(xcms)
library(RANN)
library(Biobase)
library(seqinr)
library(plyr)
library(lubridate)
library(timeDate)
library(ggplot2)
library(stringr)
require(reshape2)
require(dplyr)
library(IPO)
  
  for (j in 1:length(FractionList)){
    Fraction <- FractionList[j]
    print(Fraction)
    DataDIR <- as.character(Dirs[Fraction, "DataDIR"])
    DatFiles <- list.files(DataDIR, recursive = TRUE, full.names = TRUE)
    #DatFiles <- DatFiles[!grepl("Icon", DatFiles)]
    
    #Save your parameters used for your xcms-ing
    ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
    setwd(ResultsDIR) 
    write.csv(Params[,1:4], paste(Dirs[Fraction, "ResultsDIR"], sep = "","/Params_used.csv", collapse = NULL), row.names = TRUE)

# Getting samples and peak picking -------------------------------------------------------

setwd(DataDIR)

xset <- xcmsSet(
  DatFiles, 
  method = "centWave",  
  ppm= Params["PPM", Fraction], 
  peakwidth= c(Params["PEAKWIDTHlow", Fraction], Params["PEAKWIDTHhigh", Fraction]),
  snthresh= Params["SNthresh", Fraction], 
  mzCenterFun="apex", 
  prefilter= c(Params["PREFILTERlow", Fraction], Params["PREFILTERhigh", Fraction]),
  integrate = 1, 
  fitgauss= TRUE)


print("peak picking complete")

# Initial peak alignment -------------------------------------------------------

xset <- group(
  xset, 
  method="density", 
  bw= Params["BW", Fraction], 
  minfrac= Params["MINFRAC", Fraction],
  minsamp= Params["MINSAMP", Fraction], 
  mzwid= Params["MZWID", Fraction], 
  max= Params["MAX", Fraction])

setwd(ResultsDIR)
save(xset, file=paste(Fraction, "xset.RData", sep="."))

print(paste(Fraction,"Initial peak alignment complete", sep=" "))

# RT correction and grouping (of some iterations) -------
for (k in 1:Params["RTITs", Fraction]){

xset2 <- retcor(
  xset, 
  method = "peakgroups",
  missing = Params["MISSING", Fraction],
  extra = Params["EXTRA", Fraction], 
  smooth = "loess",
  family = "symmetric", 
  plottype = NULL
)


print(paste (Fraction, "RT correction complete", sep=" "))

# Peak align the RT-corrected data

xset2 <- group(
  xset2, 
  method="density", 
  minfrac= Params["MINFRAC", Fraction],
  minsamp= Params["MINSAMP", Fraction],
  mzwid= Params["MZWID", Fraction], 
  max= Params["MAX", Fraction],
  bw= Params["BW", Fraction])

}

xset.unfilled <- peakTable(xset2)

setwd(ResultsDIR)
write.csv(xset.unfilled, file=paste(Fraction, "Unfilled_peaks.csv", sep="."))
save(xset2, file=paste(Fraction, "xset2.RData", sep="."))

print(paste(Fraction, "additional peak alignment done", sep=""))

# Recursive peak filling -------------------------------------------------------

xset3 <- fillPeaks(xset)

setwd(ResultsDIR)
save(xset3, file=paste(Fraction, "xset3.RData", sep="."))

#xset3 is the input file for CAMERA or diffreport. 

print(paste(Fraction, "recursive peak filling is done", sep=" "))

#Calculate Peak Picking Score----

setwd(ResultsDIR)
PPS <- as.data.frame(calcPPS(xset3))
write.csv(PPS, "PPS_results.csv")

#RT Correction plot----
png("RTCorrplot.png", width = 4, height = 5, units = "in", res = 300) 
plotrt(xset2, leg = F, densplit = T)
dev.off() 
save(xset2,  file="xset2_RTCorrected.RData")

# Generate a data.frame with all the peaks----

xset.allpeaks <- peakTable(xset3)
xset.allpeaks$MassFeature <- paste("I", round((xset.allpeaks$mz),digits=4), 
                                       "R", round( xset.allpeaks$rt/60, digits=2), sep="")
xset.allpeaks$groupname <- groupnames(xset3)
xset.allpeaks$RT <- xset.allpeaks$rt/60

setwd(ResultsDIR)
write.csv(xset.allpeaks, file=paste(Fraction, "Allpeaks.table","csv",sep="."))

print(paste(Fraction, "allpeaks data table is saved", sep=" "))

#XCMSlist <- list(xset3, xset.allpeaks)
#return(invisible(XCMSlist))
}
}



