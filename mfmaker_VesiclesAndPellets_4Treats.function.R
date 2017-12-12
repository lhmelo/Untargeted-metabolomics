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

mfmaker_ttest <- function(pvalue){
  
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
  
  Blk2 <- names[grepl("_Blk_2", names)]
  Blk1 <- names[grepl("_Blk_1", names)]
  ProcessBlk <- names[grepl("ProcessBlk", names)]
  
  # Treat1ID <- "Pro9313"
  # Treat2ID <- "Pro9312"
  # Treat3ID <- "Vesicle9313"
  # Treat4ID <- "Vesicle9312"
  
  Pellets9313 <- Samples[grepl(as.character(Treat1ID), Samples)]
  Pellets9312 <- Samples[grepl(as.character(Treat2ID), Samples)]
  Vesicles9313 <- Samples[grepl(as.character(Treat3ID), Samples)]
  Vesicles9312 <- Samples[grepl(as.character(Treat4ID), Samples)]
  
  #calculate averages for each group
  # xset.allpeaks$AvePellets <- rowMeans(xset.allpeaks[, Pellets])
  # xset.allpeaks$AveVesicles <- rowMeans(xset.allpeaks[, Vesicles])
  # xset.allpeaks$AveBlk1 <- rowMeans(xset.allpeaks[, Blk1])
  # xset.allpeaks$AveBlk2 <- rowMeans(xset.allpeaks[, Blk2])
  # xset.allpeaks$AveProcessBlk <- rowMeans(xset.allpeaks[, ProcessBlk])
  
  xset.allpeaks <- xset.allpeaks %>% 
    mutate(AvePro9313 = rowMeans(xset.allpeaks[, Pellets9313])) %>%
    mutate(AvePro9312 = rowMeans(xset.allpeaks[, Pellets9312])) %>%
    mutate(AveVesicles9313 = rowMeans(xset.allpeaks[, Vesicles9313])) %>%
    mutate(AveVesicles9312 = rowMeans(xset.allpeaks[, Vesicles9312])) %>%
    mutate(AveBlk1 = rowMeans(xset.allpeaks[, Blk1])) %>%
    mutate(AveBlk2 = rowMeans(xset.allpeaks[, Blk2])) %>%
    mutate(AveProcessBlk = rowMeans(xset.allpeaks[, ProcessBlk])) %>%
    mutate(Pro9313.sd = apply(xset.allpeaks[, Pellets9313], 1, sd)) %>%
    mutate(Pro9312.sd = apply(xset.allpeaks[, Pellets9312], 1, sd)) %>%
    mutate(Vesicles9313.sd = apply(xset.allpeaks[, Vesicles9313], 1, sd)) %>%
    mutate(Vesicles9312.sd = apply(xset.allpeaks[, Vesicles9312], 1, sd)) %>%
    mutate(Blk1.sd = apply(xset.allpeaks[, Blk1], 1, sd)) %>%
    mutate(Blk2.sd = apply(xset.allpeaks[, Blk2], 1, sd)) %>%
    mutate(ProcessBlk.sd = apply(xset.allpeaks[, ProcessBlk], 1, sd))

  
  #Specify treatment groups within samples

 
  #setwd(outputpath)
  #save(xset.allpeaks, file=paste(extractiontype, "xset.allpeaks.RData", sep="."))
  #write.csv(xset.allpeaks, file=paste(extractiontype, "Allpeaks.table","csv",sep="."))
  
  
  # Retaining only mass features that elute between RTMIN and RTMAX, defined in Params
  xset.filtered <- subset(xset.allpeaks, 
                          xset.allpeaks$rt > Params["RTMIN", Fraction] & xset.allpeaks$rt <                           Params["RTMAX", Fraction])
  
  #Kepp peaks that are at least twice as big as both respective blanks in one average sample and for which all other samples are within one sd of the average blanks
  
  
  
  # xset.filtered <- subset(
  #   xset.filtered,
  #   ((xset.filtered$AveBlk1 * 2) < xset.filtered$AveVesicles &
  #      (xset.filtered$AveBlk2 * 2) < xset.filtered$AveVesicles &
  #     (xset.filtered$AveProcessBlk - xset.filtered$ProcessBlk.sd) < (xset.filtered$AvePellet + xset.filtered$Pellets.sd) &
  #        (xset.filtered$AveBlk2  - xset.filtered$Blk2.sd) < (xset.filtered$AvePellet + xset.filtered$Pellets.sd)) | 
  # 
  # ((xset.filtered$AveProcessBlk * 2) < xset.filtered$AvePellet &
  #   (xset.filtered$AveBlk2 * 2) < xset.filtered$AvePellet &
  #   (xset.filtered$AveBlk1 - xset.filtered$Blk1.sd) < (xset.filtered$AveVesicles + xset.filtered$Vesicles.sd) &
  #   (xset.filtered$AveBlk2  - xset.filtered$Blk2.sd) < (xset.filtered$AveVesicles + xset.filtered$Vesicles.sd))
  # 
  # )
  xset.filtered_test <- xset.filtered %>%
    filter((AveVesicles9313 > (AveBlk1 * 2) & 
             AveVesicles9313 > (AveBlk2 * 2) & 
             (AveVesicles9312 + Vesicles9312.sd) > (AveBlk1 - Blk1.sd) &
             (AveVesicles9312 + Vesicles9312.sd) > (AveBlk2 - Blk2.sd) &
             (AvePro9312 + Pro9312.sd) > (AveBlk2 - Blk2.sd) &
             (AvePro9312 + Pro9312.sd) > (AveProcessBlk - ProcessBlk.sd) &
             (AvePro9313 + Pro9313.sd) > (AveBlk2 - Blk2.sd) &
             (AvePro9313 + Pro9313.sd) > (AveProcessBlk - ProcessBlk.sd)) |
             
             (AveVesicles9312 > (AveBlk1 * 2) & 
                AveVesicles9312 > (AveBlk2 * 2) & 
                (AveVesicles9313 + Vesicles9313.sd) > (AveBlk1 - Blk1.sd) &
                (AveVesicles9313 + Vesicles9313.sd) > (AveBlk2 - Blk2.sd) &
                (AvePro9312 + Pro9312.sd) > (AveBlk2 - Blk2.sd) &
                (AvePro9312 + Pro9312.sd) > (AveProcessBlk - ProcessBlk.sd) &
                (AvePro9313 + Pro9313.sd) > (AveBlk2 - Blk2.sd) &
                (AvePro9313 + Pro9313.sd) > (AveProcessBlk - ProcessBlk.sd)) |
  
             
             (AvePro9313 > (AveProcessBlk * 2) & 
                AvePro9313 > (AveBlk2 * 2) & 
                (AveVesicles9312 + Vesicles9312.sd) > (AveBlk1 - Blk1.sd) &
                (AveVesicles9312 + Vesicles9312.sd) > (AveBlk2 - Blk2.sd) &
                (AveVesicles9313 + Vesicles9313.sd) > (AveBlk1 - Blk1.sd) &
                (AveVesicles9313 + Vesicles9313.sd) > (AveBlk2 - Blk2.sd) &
                (AvePro9312 + Pro9312.sd) > (AveBlk2 - Blk2.sd) &
                (AvePro9312 + Pro9312.sd) > (AveProcessBlk - ProcessBlk.sd)) |
             
             (AvePro9312 > (AveProcessBlk * 2) & 
                AvePro9312 > (AveBlk2 * 2) & 
                (AveVesicles9312 + Vesicles9312.sd) > (AveBlk1 - Blk1.sd) &
                (AveVesicles9312 + Vesicles9312.sd) > (AveBlk2 - Blk2.sd) &
                (AveVesicles9313 + Vesicles9313.sd) > (AveBlk1 - Blk1.sd) &
                (AveVesicles9313 + Vesicles9313.sd) > (AveBlk2 - Blk2.sd) &
                (AvePro9313 + Pro9313.sd) > (AveBlk2 - Blk2.sd) &
                (AvePro9313 + Pro9313.sd) > (AveProcessBlk - ProcessBlk.sd)) 
           )
             
             
             
             
  Treatsdf <- xset.filtered[, c(Pellets9313, Vesicles9313) ]
  
  #add a Pvalue for between the two treatments for QC
  xset.filtered$pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[1:3],x[4:6])$p.value}) 
  
  
  setwd(ResultsDIR)
  write.csv(xset.filtered, file=paste(Fraction, "xset.filtered.csv", sep="."), row.names = FALSE)
  
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
  write.csv(MFs, paste(Fraction, "MFs.csv", sep="."), row.names = FALSE)
  
  #List----
  
  #MFList <- list(mfs=MFs, sig=sig.groupnames, b50=bottom.fifty, filtered=xset.filtered)
  #return(invisible(MFList))
  }  
}