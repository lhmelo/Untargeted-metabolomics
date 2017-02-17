summarystatwriter <- function() {
  
  for (j in 1:length(FractionList)) {
    Fraction <- FractionList[j]
    print(Fraction)
    
    ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
    setwd(ResultsDIR)
    xset.allpeaks <-
      read.csv(paste(Fraction, "Allpeaks.table", "csv", sep = "."))
    SummaryStats["MFs_allpeaks", Fraction] <- length(xset.allpeaks$mz)
    
    xset.filtered <- read.csv(paste(Fraction, "xset.filtered.csv", sep="."))
    SummaryStats["MFs_filtered", Fraction] <- length(xset.filtered$mz)
    
    NeutralMassList <- read.csv(paste(Fraction, "NeutralMassList.csv", sep="."))
    SummaryStats["MFs_NML", Fraction] <- length(NeutralMassList$mz)
    
    UniqueNML <- read.csv(paste(Fraction, "UniqueNML.csv", sep="."))
    SummaryStats["MFs_UniqueNML", Fraction] <- length(UniqueNML[,1])
    
    UniqueHypNML <- read.csv(paste(Fraction, "UniqueHypNM.csv", sep="."))
    SummaryStats["MFs_UniqueHypNML", Fraction] <- length(UniqueHypNML[,1])
    
    load("sig.groupnames.RData")
    SummaryStats["sig.peaks", Fraction] <- length(sig.groupnames[,1])  
  }
  
    SummaryStats$percentAq <- (SummaryStats$CyanoAq/(SummaryStats$CyanoAq + SummaryStats$CyanoDCM))*100
    
    SummaryStats$percentDCM <- (SummaryStats$CyanoDCM/(SummaryStats$CyanoAq + SummaryStats$CyanoDCM))*100
    
    setwd(ExpDIR)
    write.csv(SummaryStats, file="SummaryStats.csv")
    
  return(SummaryStats)
}