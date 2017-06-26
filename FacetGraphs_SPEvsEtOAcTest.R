
require(stringr)
require(plyr)
require(cowplot)
require(dplyr)
require(ggplus)

#data <- VibrioDCM
#class(data$experiment)
#data$experiment <-as.factor(data$experiment)
#class(data$experiment)

#description <-"AQ"
#column <- "Cyano"
#culture <- "Vibrio"

FacetGraphs <- function(data, column, description) {
  
  
    data.summ <- data 
    
    p <- ggplot(data.summ, aes(y=mean.Area, x=SampID)) + 
      geom_bar(stat="identity", fill="white", colour="black", width = 0.5, position = "dodge") +
      geom_errorbar(aes(ymin=mean.Area-sd.Area, ymax=mean.Area+sd.Area), width=0.2) +
      #labs(y="Area", x="SampID") +
      scale_y_continuous(expand = c(0,0))+
      theme(axis.text.x=element_text(angle=50, size=8, vjust=0.5)) +
      facet_wrap(~Compound.Name, scales= "free", ncol=3)
    
    #facet_multiple(plot = p, 
    #               facets = 'Compound.Name', 
    #               ncol = 2, 
    #               nrow = 2,
    #               scales= "free")
    
    return(p)
    
  }

  
  
  
  
  
