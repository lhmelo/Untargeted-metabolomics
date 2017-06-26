
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
  
  
    data.summ <- data %>%
      group_by(SampID, Compound.Name) %>%
      dplyr::summarise(
        mean.Area = mean(Area, na.rm = TRUE),
        sd.Area = sd(Area, na.rm = TRUE))
    
    #data.summ.trim <- data.summ[!is.na(data.summ$mean.Area),]
    
    p <- ggplot(data.summ %>% filter(SampID == "Vesicle9313" | SampID == "Vesicle9312"), aes(y=mean.Area, x=SampID)) + 
      geom_bar(stat="identity", fill="white", colour="black", width = 0.5, position = "dodge") +
      geom_errorbar(aes(ymin=mean.Area-sd.Area, ymax=mean.Area+sd.Area), width=0.2) +
      labs(y="Area", x="SampID") +
      scale_y_continuous(expand = c(0,0))+
      facet_wrap(~Compound.Name, scales= "free", ncol=3 )
    
    #facet_multiple(plot = p, 
    #               facets = 'Compound.Name', 
    #               ncol = 2, 
    #               nrow = 2,
    #               scales= "free")
    
    return(p)
    
  }

  
  
  
  
  
