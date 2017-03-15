#logplot()  
#Input:
#1. xset.filtered, which has already been filtered against the blanks.
#
#Output:
#1. logplotdata.csv: data derived from xset.filtered which has been filtered to remove rows where p > 0.05.
#
#logplot() performs a t test for each mass feature, retains only mass features which display signficant distributions between samples (p < 0.05) and calculates a new "logratio".  


logplot <- function (sigcutoff = 0.05){
  
  for (j in 1:length(FractionList)) {
   #data <- NA
   #pvalue <- NA
   #names <- NA
   #Samples <- NA
   #Treat1 <- NA
   #Treat2 <- NA
   #sig.groupnames <- NA
   
    Fraction <- FractionList[j]
    print(Fraction)
    
    ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
    setwd(ResultsDIR)
    data <-
      read.csv(paste(Fraction, "xset.filtered.csv", sep = "."))
    
    
    
    #Treat1ID <- "9313"
    #Treat2ID <- "9312"
    
    names <- colnames(data)
    Samples <- names[grepl("Smp", names)]
    #Blanks <- names[grepl("Blk", names)]
    #Pooled <- names[grepl("Poo", names)]
    
    Treat1 <- Samples[grepl(as.character(Treat1ID), Samples)]
    Treat2 <- Samples[grepl(as.character(Treat2ID), Samples)]
    
    #Which are the most abundant features?
    #Are any of these features in common between the two strains?
    
    data$AveTreat1 <- rowMeans(data[, Treat1])
    data$AveTreat2 <- rowMeans(data[, Treat2])
    #test <- data[with(data, order(-data$Ave9313)),]
    
    topTreat1 <- data %>% arrange(desc(AveTreat1)) %>% head(50)
    #top9313 <- top_n(data, 50, data$Ave9313)
    write.csv(topTreat1,
              file = paste(Fraction, "top", Treat1ID, "csv", sep = "."),
              row.names = FALSE)
    #          
    topTreat2 <- data %>% arrange(desc(AveTreat2)) %>% head(50)         
    #top9312 <- top_n(data, 50, data$Ave9312)
    write.csv(topTreat2,
              file = paste(Fraction, "top", Treat2ID, "csv", sep = "."),
              row.names = FALSE)
    
    #Are any of these features in common between the two strains?
    #top9313_9312 <- intersect(top9312$groupname, top9313$groupname)
    #
    
   
    
    #building "logratio" which will represent enrichment of 9313 mass features as positive values and encrichment of 9312 mass features as negative values
    #
    #sig.groupnames$logratio <- NA
    
    data$log2ratio <-
      ifelse(
       data$AveTreat1 > data$AveTreat2,
        log2(data$AveTreat1 / data$AveTreat2),
        ifelse(
          data$AveTreat1 < data$AveTreat2,
          -(log2(
            data$AveTreat2 / data$AveTreat1
          )),
          "other"
        )
      )
    
   data$log2ratio <- gsub('Inf', '25', data$log2ratio)
    data$log2ratio <- as.numeric(data$log2ratio)
    
    sig.groupnames <- subset(data, data$pvalue < sigcutoff)
    # sig.groupnames$AveTreat1 <- rowMeans(sig.groupnames[, Treat1])
    # sig.groupnames$AveTreat2 <- rowMeans(sig.groupnames[, Treat2])
    
    #replacing Inf with "25" in order to treat the column as numeric
    sig.groupnames$log2ratio <- gsub('Inf', '25', sig.groupnames$log2ratio)
    sig.groupnames$log2ratio <- as.numeric(sig.groupnames$log2ratio)
    
    write.csv(data, paste(Fraction, "logplotdata", "csv", sep="."))
    
    p <- ggplot(sig.groupnames, aes(x = mz, y = log2ratio, color = pvalue)) +
      geom_point() +
      scale_color_continuous(low = "red", high = "blue") +
      labs(title = (paste(Fraction, Experiment, "Raw Data", sep = " "))) +
      theme_classic()
    
    print(p)
    
    save_plot(paste(Fraction, Experiment, "RawDataLogPlot.pdf", sep="."), p,  base_aspect_ratio = 2)
}
}