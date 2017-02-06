# mfmaker_ANOVA (Vibrio Samples Edition)----
# 
# This function is tailored to an experiment compairing metabalomes of vibrio fischeri grown under three different conditions and uses ANOVA to define lists of significant peaks to be used in QC plots and refining neutral mass lists
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


mfmaker_ANOVA <- function(pvalue){
  
  #Basic filtering-----
  
  for (j in 1:length(FractionList)){
    Fraction <- FractionList[j]
    print(Fraction)
    ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
    setwd(ResultsDIR)
    xset.allpeaks <- read.csv(paste(Fraction, "Allpeaks.table","csv",sep="."))
    
    names <- colnames(xset.allpeaks)
 #Specify blanks, samples and pooled samples.
  Samples <-names[grepl("Smp", names)]
  Blanks <- names[grepl("Blk", names)]
  Pooled <- names[grepl("Poo", names)]
  
  #calculate averages for each group
  xset.allpeaks$AveSmp<- rowMeans(xset.allpeaks[, Samples])
  xset.allpeaks$AveBlk<- rowMeans(xset.allpeaks[, Blanks])
  xset.allpeaks$AvePoo<- rowMeans(xset.allpeaks[, Pooled])
  
  #Specify treatment groups within samples
  
  Treat1 <- Samples[grepl(as.character(Treat1ID), Samples)]
  Treat2 <- Samples[grepl(as.character(Treat2ID), Samples)]
    
 
  #setwd(outputpath)
  #save(xset.allpeaks, file=paste(extractiontype, "xset.allpeaks.RData", sep="."))
  #write.csv(xset.allpeaks, file=paste(extractiontype, "Allpeaks.table","csv",sep="."))
  
  
  # Retaining only mass features that elute between RTMIN and RTMAX, defined in Params
  xset.filtered <- subset(xset.allpeaks, 
                          xset.allpeaks$rt > Params["RTMIN", Fraction] & xset.allpeaks$rt <                           Params["RTMAX", Fraction])
  
  #Get rid of any peaks that are bigger in the blank than in the average pooled sample
  
  xset.filtered <- subset(xset.filtered, xset.filtered$AveBlk < xset.filtered$AvePoo) 
  Treatsdf <- xset.filtered[, c(Treat1,Treat2) ]
  
  #add a Pvalue for between the two treatments for QC
  xset.filtered$pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[1:3],x[4:6])$p.value}) 
  
  
  setwd(ResultsDIR)
  write.csv(xset.filtered, file=paste(Fraction, "xset.filtered.csv", sep="."))
  
  #Pull out top 50 most differentially expressed peaks

  sig.groupnames <- subset(xset.filtered, xset.filtered$pvalue < 0.05)
  
  save(sig.groupnames, file="sig.groupnames.RData")
  
  bottom.fifty <- top_n(sig.groupnames,-50, pvalue)
  
  #Quality control ---- 
  
  #Aq.Samples <- Aq.Samples[!grepl("Poo|Blk", Aq.Samples)]
  #Aq7 <- Aq.Samples[grepl("7", Aq.Samples)]
  #Aq10 <- Aq.Samples[grepl("10", Aq.Samples)]
  #AQ12 <- Aq.Samples[grepl("12", Aq.Samples)]
  
  #next line: this is where you could filter samples according to p-values generated in xset.filtered (in place of runif, which selects random samples)
  #MFs <- xset.filtered$groupname[runif(30, min = 1, max = length(xset.filtered$groupname))] 
  MFs <- bottom.fifty$groupname
  MFs<- as.character(MFs)
  setwd(ResultsDIR)
  write.csv(MFs, paste(Fraction, "MFs.csv", sep="."))
  
  #List----
  
  #MFList <- list(mfs=MFs, sig=sig.groupnames, b50=bottom.fifty, filtered=xset.filtered)
  #return(invisible(MFList))
  }  
}