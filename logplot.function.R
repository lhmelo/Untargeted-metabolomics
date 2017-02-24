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
    
    topTreat1 <- data %>% arrange(desc(data$AveTreat1)) %>% head(50)
    #top9313 <- top_n(data, 50, data$Ave9313)
    write.csv(topTreat1,
              file = paste(Fraction, "top", Treat1ID, "csv", sep = "."),
              row.names = FALSE)
    #          
    topTreat2 <- data %>% arrange(desc(data$AveTreat2)) %>% head(50)         
    #top9312 <- top_n(data, 50, data$Ave9312)
    write.csv(topTreat2,
              file = paste(Fraction, "top", Treat2ID, "csv", sep = "."),
              row.names = FALSE)
    
    #Are any of these features in common between the two strains?
    #top9313_9312 <- intersect(top9312$groupname, top9313$groupname)
    #
    
    sig.groupnames <- subset(data, data$pvalue < sigcutoff)
    sig.groupnames$AveTreat1 <- rowMeans(sig.groupnames[, Treat1])
    sig.groupnames$AveTreat2 <- rowMeans(sig.groupnames[, Treat2])
    
    #building "logratio" which will represent enrichment of 9313 mass features as positive values and encrichment of 9312 mass features as negative values
    #
    #sig.groupnames$logratio <- NA
    
    sig.groupnames$logratio <-
      ifelse(
        sig.groupnames$AveTreat1 > sig.groupnames$AveTreat2,
        log10(sig.groupnames$AveTreat1 / sig.groupnames$AveTreat2),
        ifelse(
          sig.groupnames$AveTreat1 < sig.groupnames$AveTreat2,
          -(log10(
            sig.groupnames$AveTreat2 / sig.groupnames$AveTreat1
          )),
          "other"
        )
      )
    
    #replacing Inf with "4" in order to treat the column as numeric
    sig.groupnames$logratio <- gsub('Inf', '5', sig.groupnames$logratio)
    sig.groupnames$logratio <- as.numeric(sig.groupnames$logratio)
    
    
    p <- ggplot(sig.groupnames, aes(x = mz, y = logratio, color = pvalue)) +
      geom_point() +
      scale_color_continuous(low = "red", high = "blue") +
      labs(title = (paste(Fraction, Experiment, sep = " "))) +
      theme_classic()
    
    print(p)
}
}