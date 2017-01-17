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