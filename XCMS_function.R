
# xcmsfunction----
# This script uses XCMS to generate a list of peaks that have been
# aligned, RT corrected, realigned, and recursively filled. 
# 
# Before running the function:
# 1. Define outputpath and mzxmlpath
# 2. Define extractiontype in function call
# 
# Input:
# 1. .mzXML sample files unbinned in a single folder
# 
# Output is a list containing:
# 1. a processed xcmsSet object (xset3)
# 2. a data.frame containing a list of peaks that have been
# aligned, RT corrected, realigned, and recursively filled (xset.allpeaks)

xcmsfunction <- function(extractiontype) {

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
  


# Getting samples and peak picking -------------------------------------------------------

setwd(mzxmlpath)

Samples <- list.files( getwd(), pattern=".mzXML", full.names=F, 
                      recursive = TRUE) # looking for appropriate files

sampleset <- Samples[grepl(extractiontype, Samples)]
print(sampleset)

xset <- xcmsSet(
  sampleset, 
  method = "centWave",  
  ppm= PPM, 
  peakwidth= PeakWidth,
  snthresh= SNthresh, 
  mzCenterFun="apex", 
  prefilter= Prefilter ,
  integrate = 1, 
  fitgauss= TRUE)


print("peak picking complete")

# Initial peak alignment -------------------------------------------------------

xset <- group(
  xset, 
  method="density", 
  bw= BW, 
  minfrac= MinFrac,
  minsamp= Minsamp, 
  mzwid= MzWid, 
  max= MAX)

setwd(outputpath)
save(xset, file=paste(extractiontype, "xset.RData", sep="."))

print("Initial peak alignment complete")

# RT correction -------------------------------------------------------


xset2 <- retcor(
  xset, 
  method = "peakgroups",
  missing = Missing, 
  extra = Extra, 
  smooth = "loess",
  family = "symmetric", 
  plottype = NULL
)


print("RT correction complete")

# Peak align the RT-corrected data -------------------------------------------------------

xset2 <- group(
  xset2, 
  method="density", 
  minsamp=Minsamp, 
  minfrac=MinFrac, 
  mzwid= MzWid, 
  bw= BW, 
  max=MAX)


xset.unfilled <- peakTable(xset)

setwd(outputpath)
write.csv(xset.unfilled, file=paste(extractiontype, "Unfilled_peaks.csv", sep="."))
save(xset2, file=paste(extractiontype, "xset2.RData", sep="."))

print("second peak alignment done")

# Recursive peak filling -------------------------------------------------------

xset3 <- fillPeaks(xset)

setwd(outputpath)
save(xset3, file=paste(extractiontype, "xset3.RData", sep="."))

print("recursive peak filling is done")

#xset3 is the input file for CAMERA or diffreport. 

# Generate a data.frame with all the peaks----

xset.allpeaks <- peakTable(xset3)
xset.allpeaks$MassFeature <- paste("I", round((xset.allpeaks$mz),digits=4), 
                                       "R", round( xset.allpeaks$rt/60, digits=2), sep="")
xset.allpeaks$groupname <- groupnames(xset3)
xset.allpeaks$RT <- xset.allpeaks$rt/60

setwd(outputpath)
write.csv(xset.allpeaks, file=paste(extractiontype, "Allpeaks.table","csv",sep="."))

print("allpeaks data table is saved")

XCMSlist <- list(xset3, xset.allpeaks)
return(invisible(XCMSlist))

}


# mfmaker_ANOVA----
# 
# This function is tailored to Eddy Cruise data and uses ANOVA to define lists of significant peaks to be used in QC plots and refining neutral mass lists
# 
# Input is 
# 1. xset.allpeaks
# 2. extractiontype
# 3. A pvalue used to set a limit for filtering ANOVA output (pvalue)
# 
# Output is
# 1. a list of Mass Features (MFs) to serve as input for QC plots
# 2. a subset of allpeaks which ANOVA determines are significant at a level determined by  pvalue (sig.groupnames)
# 3. the fifty most significant peaks from sig.groupnames, to be plotted in QC step (bottom.fifty)
# 4. xset.filtered

mfmaker_ANOVA <- function(xset.allpeaks, extractiontype, pvalue){

#Basic filtering-----

#xset.allpeaks$AveBlank<- (xset.allpeaks$X160806_Blk_8.3_FilterBlk + 
#                                xset.allpeaks$X160806_Smp_BAB)/2
xset.allpeaks$AveSmp<- rowMeans(xset.allpeaks[, 9:17])

setwd(outputpath)
save(xset.allpeaks, file=paste(extractiontype, "xset.allpeaks.RData", sep="."))
write.csv(xset.allpeaks, file=paste(extractiontype, "Allpeaks.table","csv",sep="."))


# Retaining only mass features that elute after 0.5 min
xset.filtered <- subset(xset.allpeaks, xset.allpeaks$rt > 30)

if(extractiontype == "DCM") {
  xset.filtered <- subset(xset.filtered, xset.filtered$FilterBlk_DCM < xset.filtered$AveSmp)
} else if (extractiontype == "Aq") {
xset.filtered <- subset(xset.filtered, xset.filtered$FilterBlk_Aq < xset.filtered$AveSmp)
} else {
  stop()
}

setwd(outputpath)
write.csv(xset.filtered, file=paste(extractiontype, "xset.filtered", sep="."))

if(extractiontype == "DCM")
{
  xset.filtered.trim <- xset.filtered %>%
    select(X7A_DCM,
           X7B_DCM,
           X7C_DCM,
           X10A_DCM,
           X10B_DCM,
           X10C_DCM,
           X12A_DCM,
           X12B_DCM,
           X12C_DCM,
           groupname)
} else if(extractiontype == "Aq")
{
  xset.filtered.trim <- xset.filtered %>%
    select(X7A_Aq,
           X7B_Aq,
           X7C_Aq,
           X10A_Aq,
           X10B_Aq,
           X10C_Aq,
           X12A_Aq,
           X12B_Aq,
           X12C_Aq,
           groupname)
} else {
 stop(); 
}



#Transpose and melt data.frame----
transposed.df <- t(xset.filtered.trim)
colnames(transposed.df) <- transposed.df["groupname",]
transposed.df <- transposed.df[-10,]
transposed.df <- data.frame(transposed.df)

setwd(outputpath)
write.csv(transposed.df, file = paste(extractiontype, "transposed.df.csv", sep="."))

transposed.df$Samples <- rownames(transposed.df)
transposed.df$Stations <- c(7,7,7,10,10,10,12,12,12)

#move Sample names and Stations to front of data.frame
D <- dim(transposed.df)
transposed.df <- transposed.df[,c(length(transposed.df),1:length(transposed.df)-1)]
transposed.df <- transposed.df[,c(length(transposed.df),1:length(transposed.df)-1)]

melt.df <-
  melt(
    transposed.df, id.vars = c("Samples","Stations"), variable.name = "groupname", value.name =
      "Area"
  )

#ANOVA----

ANOVAfunction <- function(data)
{
  aov(Area ~ Stations, data = data)
  
}


anova.output <- dlply (melt.df, .(groupname), ANOVAfunction)

pvaluepicker <- function(x) {
  c(summary(x)[[1]][["F value"]][[1]],
    summary(x)[[1]][["Pr(>F)"]][[1]])
}


anova.summary <- ldply(anova.output, pvaluepicker)
anova.summary <- data.frame(anova.summary)

colnames(anova.summary)[2] <- 'Fvalue'
colnames(anova.summary)[3] <- 'Pvalue'

sig.groupnames <-
  subset(anova.summary, anova.summary$Pvalue < pvalue)

save(sig.groupnames, file="sig.groupnames.RData")

bottom.fifty <- top_n(sig.groupnames,-50, Pvalue)

#Quality control ---- 

#Aq.Samples <- Aq.Samples[!grepl("Poo|Blk", Aq.Samples)]
#Aq7 <- Aq.Samples[grepl("7", Aq.Samples)]
#Aq10 <- Aq.Samples[grepl("10", Aq.Samples)]
#AQ12 <- Aq.Samples[grepl("12", Aq.Samples)]

#next line: this is where you could filter samples according to p-values generated in xset.filtered (in place of runif, which selects random samples)
#MFs <- xset.filtered$groupname[runif(30, min = 1, max = length(xset.filtered$groupname))] 
MFs <- bottom.fifty$groupname
MFs<- as.character(MFs)
setwd(outputpath)
write.csv(MFs, paste(extractiontype, "xset_50mostsignif_MassFeatures.csv", sep="."))

#List----

MFList <- list(mfs=MFs, sig=sig.groupnames, b50=bottom.fifty, filtered=xset.filtered)
return(invisible(MFList))

}


#samplesetmaker----
#This function defines a group of samples based on the extraction phase. This function is used to define the sampleset variable used in the subsequent function, qcplotter.

samplesetmaker <- function(extractiontype){
  
  setwd(mzxmlpath)
  
  Samples <- list.files( getwd(), pattern=".mzXML", full.names=F, 
                         recursive = TRUE) # looking for appropriate files
  
  sampleset <- Samples[grepl(extractiontype, Samples)]
  
  return(sampleset)
}

qcplotter <- function(xset, sampleset, xset.filtered, MFs){

EIC.uncorrected <- list()
EIC.corrected <- list()

# This next step will take some time to process, so don't expect instant results. 
setwd(mzxmlpath)
RandSamp <- round(runif(1, min = 1, max = length(sampleset)))

#for (i in 1:length(MFs)){   
#  EIC.uncorrected[[i]] <- getEIC(xset, rt="raw", groupidx=MFs[i], sampleidx=Aq.Samples)
#  EIC.corrected[[i]] <- getEIC(xset, rt="corrected", groupidx=MFs[i], sampleidx=Aq.Samples)
#}

for (i in 1:length(MFs)){   
  EIC.uncorrected[[i]] <- getEIC(xset, rt="raw", groupidx=MFs[i])
  EIC.corrected[[i]] <- getEIC(xset, rt="corrected", groupidx=MFs[i])
}



ColRainbow <- colorRampPalette(c("green", "blue", "purple")) 
MyColors <- c(ColRainbow(RandSamp-1), "red")

#setwd(mzxmlpath)

#
LastSamp <- sampleset[RandSamp[length(RandSamp)]]

xset.raw <- xcmsRaw(LastSamp, profstep=0.01, profmethod="bin")


setwd(outputpath)
pdf(paste(Sys.Date(), extractiontype, "xset_QC.pdf", sep = "."), 8.5,11)

par(mfrow=c(4,3), mar=c(3,3,3,0.5)) 
for(i in 1:length(MFs)){
  plot(EIC.uncorrected[[i]], xset, groupidx=1, rtrange=60, col=MyColors, main=MFs[i])
  mtext(paste(i, xset.filtered$MassFeature[xset.filtered$groupname == MFs[i]]), 
        side=3, line=-1, adj=0, padj=0, cex=0.8)
  plot(EIC.corrected[[i]], xset, groupidx=1, rtrange=60, col=MyColors)
  RT <- xset.filtered$rt[xset.filtered$groupname == MFs[i]] 
  RTRange <- c(RT-30, RT+30)
  mz <- xset.filtered$mz[xset.filtered$groupname == MFs[i]] 
  mzRange <- c(mz-0.02, mz+0.02)
  mzRange.poly.low <- mz- mz*(0.5*PPM)/1e6
  mzRange.poly.up <- mz*(0.5*PPM)/1e6 + mz
  plotRaw(xset.raw, mzrange=mzRange, rtrange=RTRange, log=FALSE) 
  abline(h=mz, lty=2, col="gray35")
  #mtext(paste("abund =", round(xset.filtered[
  # xset.filtered$groupname == MFs[i]])), 
  #SampCol[RandSamp[length(RandSamp)]]], digits=0)), 
  #   side=3, line=-1, adj=0, padj=0, cex=0.8)
  polygon(c(RTRange[2], RTRange[1], RTRange[1], RTRange[2]), 
          c(mzRange.poly.up, mzRange.poly.up,
            mzRange.poly.low, mzRange.poly.low), 
          col=col2alpha("blue", alpha=0.1), border=NA)
  abline(v=RT, lty=2, col="gray35")
}
dev.off()

return("pdf saved")
}


