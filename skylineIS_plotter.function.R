library(cowplot)

skylineIS_plotter <- function () {
  
  for (j in 1:length(Files)){
    
    setwd(filedirectory)
    ISfile <- read.csv(Files[j], stringsAsFactors = FALSE)
    
    Fraction <- Files[j] 
    Fraction <- gsub('_AllPro.csv', '', Fraction)
    Fraction <- gsub('Ingalls_Lab_QE_Transition Results_', '', Fraction)
    
  
    data <- ISfile %>%
      filter(grepl("Vesicle", Replicate.Name)) %>%
      mutate(Replicate.Name = regmatches(Replicate.Name, regexpr("Vesicle931[1-3]_[1-3A-C]", Replicate.Name) ) ) %>%
      mutate(Area = as.numeric(Area)) %>%
      group_by(Precursor.Ion.Name)
    
   
    p <- ggplot(data, aes(y=Area, x=Replicate.Name)) + 
      geom_bar(stat="identity", fill="white", colour="black", width = 0.5, position = "dodge") +
      #geom_errorbar(aes(ymin=mean.Area-sd.Area, ymax=mean.Area+sd.Area), width=0.2) +
      labs(title = (paste(Fraction, "Skyline Standards", sep = " ")), y="Area", x="Replicate") +
      scale_y_continuous(expand = c(0,0))+
      facet_wrap(~Precursor.Ion.Name, scales= "free", ncol=3 ) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    print(p)
  
    save_plot(paste(Fraction, "SkylineISPlot.pdf", sep="."), p,  base_aspect_ratio = 2)
    
}}