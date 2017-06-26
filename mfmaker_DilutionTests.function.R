# mfmaker_ttest Vesicles and Pellets----
# 
# This function is tailored to an experiment compairing metabalomes of intercellular metabolites and extracellular vesicle metabolites and uses a t test to define lists of significant peaks to be used in QC plots and refining neutral mass lists
# 
# This function has a more rigourous blank filter than its predecesor.
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
library(tidyverse)

mfmaker_DilutionTests <- function(pvalue){
  
  #Basic filtering-----
  
  for (j in 1:length(FractionList)){
    Fraction <- FractionList[j]
    print(Fraction)
    ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
    setwd(ResultsDIR)
    xset.allpeaks <- read.csv(paste(Fraction, "Allpeaks.table","csv",sep="."))
    
    names <- colnames(xset.allpeaks)
 #Specify blanks, samples and pooled samples.
  Samples <-names[grepl("Poo", names)]
  
  Blanks <- names[grepl("_Blk_", names)]
  
  #(?!0)
  
  One_1to5 <- Samples[grepl("1mL_1to5[^0]", Samples)]
  One_1to10 <- Samples[grepl("1mL_1to10", Samples)]
  One_1to20 <- Samples[grepl("1mL_1to20", Samples)]
  One_1to50 <- Samples[grepl("1mL_1to50", Samples)]
  Hundred_FS <- Samples[grepl("100uL_FS", Samples)]
  Hundred_1to2 <- Samples[grepl("100uL_1to2", Samples)]
  Hundred_1to3 <- Samples[grepl("100uL_1to3", Samples)]
  Hundred_1to7 <- Samples[grepl("100uL_1to7", Samples)]
  
  #calculate averages for each group
  # xset.allpeaks$AvePellets <- rowMeans(xset.allpeaks[, Pellets])
  # xset.allpeaks$AveVesicles <- rowMeans(xset.allpeaks[, Vesicles])
  # xset.allpeaks$AveBlk1 <- rowMeans(xset.allpeaks[, Blk1])
  # xset.allpeaks$AveBlk2 <- rowMeans(xset.allpeaks[, Blk2])
  # xset.allpeaks$AveProcessBlk <- rowMeans(xset.allpeaks[, ProcessBlk])
  
  xset.allpeaks <- xset.allpeaks %>% 
    mutate(AveOne_1to5 = rowMeans(xset.allpeaks[, One_1to5])) %>%
    mutate(AveOne_1to10 = rowMeans(xset.allpeaks[, One_1to10])) %>%
    mutate(AveOne_1to20 = rowMeans(xset.allpeaks[, One_1to20])) %>%
    mutate(AveOne_1to50 = rowMeans(xset.allpeaks[, One_1to50])) %>%
    mutate(AveHundred_FS = rowMeans(xset.allpeaks[, Hundred_FS])) %>%
    mutate(AveHundred_1to2 = rowMeans(xset.allpeaks[, Hundred_1to2])) %>%
    mutate(AveHundred_1to3 = rowMeans(xset.allpeaks[, Hundred_1to3])) %>%
    mutate(AveHundred_1to7 = rowMeans(xset.allpeaks[, Hundred_1to7])) %>%
    mutate(AveBlanks = rowMeans(xset.allpeaks[, Blanks]))
    

  
  #Specify treatment groups within samples

 
  #setwd(outputpath)
  #save(xset.allpeaks, file=paste(extractiontype, "xset.allpeaks.RData", sep="."))
  #write.csv(xset.allpeaks, file=paste(extractiontype, "Allpeaks.table","csv",sep="."))
  
  
  # Retaining only mass features that elute between RTMIN and RTMAX, defined in Params
  xset.filtered <- subset(xset.allpeaks, 
                          xset.allpeaks$rt > Params["RTMIN", Fraction] & xset.allpeaks$rt <                           Params["RTMAX", Fraction])
  
  
  
  xset.filtered <- xset.filtered %>%
    filter(rt > Params["RTMIN", Fraction] & rt < Params["RTMAX", Fraction]) %>%
    filter((AveOne_1to5 > (AveBlanks * 2)) |
             (AveOne_1to10 > (AveBlanks * 2)) |
             (AveOne_1to20 > (AveBlanks * 2)) |
             (AveOne_1to50 > (AveBlanks * 2)) |
             (AveHundred_1to2 > (AveBlanks * 2)) |
             (AveHundred_1to3 > (AveBlanks * 2)) |
             (AveHundred_1to7 > (AveBlanks * 2)) |
             (AveHundred_FS > (AveBlanks * 2)))
  
  
  
  #Treatsdf <- xset.filtered[, c(Pellets, Vesicles) ]
  
  #add a Pvalue for between the two treatments for QC
  #xset.filtered$pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[1:3],x[4:6])$p.value}) 
  
  
  setwd(ResultsDIR)
  write.csv(xset.filtered, file=paste(Fraction, "xset.filtered.csv", sep="."), row.names = FALSE)
  
  #Pull out top 50 most differentially expressed peaks

  # sig.groupnames <- subset(xset.filtered, xset.filtered$pvalue < 0.05)
  # 
  # save(sig.groupnames, file="sig.groupnames.RData")
  # 
  # bottom.fifty <- top_n(sig.groupnames,-50, pvalue)
  
  #Quality control ---- 
  
  #Aq.Samples <- Aq.Samples[!grepl("Poo|Blk", Aq.Samples)]
  #Aq7 <- Aq.Samples[grepl("7", Aq.Samples)]
  #Aq10 <- Aq.Samples[grepl("10", Aq.Samples)]
  #AQ12 <- Aq.Samples[grepl("12", Aq.Samples)]
  
  #next line: this is where you could filter samples according to p-values generated in xset.filtered (in place of runif, which selects random samples)
  #MFs <- xset.filtered$groupname[runif(30, min = 1, max = length(xset.filtered$groupname))] 
  # MFs <- bottom.fifty$groupname
  # MFs<- as.character(MFs)
  # setwd(ResultsDIR)
  # write.csv(MFs, paste(Fraction, "MFs.csv", sep="."), row.names = FALSE)
  
  #List----
  
  #MFList <- list(mfs=MFs, sig=sig.groupnames, b50=bottom.fifty, filtered=xset.filtered)
  #return(invisible(MFList))
  }  
}