# NMLmatch function --------------------------------------------------

# This function will take two data.frames, X and Y, and output a list containing two data frames. [[1]]  contains a listing which mass features match which in the other dataset. [[2]] contains peak areas accross all samples for the mass features listed in [[1]].

# Data.frame X must include the following columns:
#       1. Compound.Name
#       2. MassFeature
#       3. mz
#       4. RT
#       
# List of target compounds should be data.frame X:
# As of 27-01-17 using "MRM_Methods_Table.csv" for Cyano and HILIC targetted compounds and labelled standards
# 
# xsAq.allpeaks object from XCMS output should be data.frame Y.
# 
# Data.frame X must include the following columns:
#       
#       1. MassFeature
#       2. mz
#       3. RT
# 
# The default is that matches must be within 5 ppm and 0.2 min, but you can
# change the defaults when you call the function, eg. with defaults:
# mfmatch(Data1, Data2)
# with other settings for PPM and RTRange:
# mfmatch(Data1, Data2, PPM = 10, RTRange = 0.3)
# column is either "RP" or "HILIC"



nmlmatch <- function(PPM = 5, RTRange = 0.2){

  require(plyr)
  require(dplyr)
  require(stringr)
  require(readr)
  
  
  for (j in 1:length(FractionList)) {
    Fraction <- FractionList[j]
    print(Fraction)
    
    ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
    setwd(ResultsDIR)
    NeutralMassList <-
      read.csv(paste(Fraction, "NeutralMassList.csv", sep="."))
    X <- NeutralMassList
    
    
    setwd(targetcompoundpath)
    data <- read.csv("TargetCompoundList.csv", stringsAsFactors = FALSE)
    setwd(ResultsDIR)
 
 POLARITY <- ifelse(Params["POLARITY", Fraction]==1, "positive",
                    ifelse(Params["POLARITY", Fraction]==2, "positive",
                           ifelse(Params["POLARITY", Fraction]==3, "positive",
                                  ifelse(Params["POLARITY", Fraction]==4, "negative", NA))))
 
 COLUMN <- ifelse(Params["POLARITY", Fraction]==1, "RP",
                    ifelse(Params["POLARITY", Fraction]==2, "RP",
                           ifelse(Params["POLARITY", Fraction]==3, "HILIC",
                                  ifelse(Params["POLARITY", Fraction]==4, "HILIC", NA))))
 

 Target.Compounds <- data[which(data[,'Column'] == COLUMN  & data[,'polarity']==POLARITY ), c("Compound.Name", "Compound.Type", "mz", "RT", "Column", "polarity")]
 
 
 Target.Compounds$mz <- as.numeric(Target.Compounds$mz)
 
 Target.Compounds <- Target.Compounds[complete.cases(Target.Compounds$mz),]
 
 Target.Compounds$MassFeature <- paste("I", round((Target.Compounds$mz),digits=4), 
                                       "R", round( Target.Compounds$RT, digits=2), sep="")
 
 
  Y <- Target.Compounds
  
  DF.X <- X[, c("MassFeature", "mz", "RT", "pcgroup", "IonType", "AvePoo", "MassOfM")]
  DF.Y <- Y[, c("Compound.Name", "MassFeature", "mz", "RT")]
  
  DF.Y$MassFeature <- as.character(DF.Y$MassFeature)
  DF.X$MassFeature <- as.character(DF.X$MassFeature)
  
  DF.X <- arrange(DF.X, MassFeature)
  DF.Y <- arrange(DF.Y, MassFeature)
  
  # Renaming to keep track of which mass feature, m/z, and RT came from 
  # which dataset.
  names(DF.X) <- paste(names(DF.X), "X", sep = ".")
  names(DF.Y) <- paste(names(DF.Y), "Y", sep = ".")
  
  MFmatch <- list()
  Matched.X <- c()
  MFname.X <- as.character(DF.X$MassFeature.X)
  mz.X <- DF.X$mz.X
  names(mz.X) <- DF.X$MassFeature.X
  RT.X <- DF.X$RT.X
  names(RT.X) <- DF.X$MassFeature.X
  
  # Checking each row in DF.X for any matches in DF.Y
  for (i in DF.X$MassFeature.X){
    
    MFmatch[[i]] <- DF.Y[
      DF.Y$mz.Y > (mz.X[i] - (PPM/1e6*mz.X[i])) 
      & DF.Y$mz.Y < (mz.X[i] + (PPM/1e6*mz.X[i])) 
      & DF.Y$RT.Y > RT.X[i] - RTRange 
      & DF.Y$RT.Y < RT.X[i] + RTRange, ]
    
    if (nrow(MFmatch[[i]]) > 0) {
      MFmatch[[i]] <- 
        rbind(MFmatch[[i]], 
              DF.Y[DF.Y$mz.Y > (mz.X[i] - 
                                  (PPM/1e6*mz.X[i])) 
                   & DF.Y$mz.Y < (mz.X[i] + 
                                    (PPM/1e6*mz.X[i])) 
                   & is.na(DF.Y$RT.Y), ])
      
      
    } else {
      MFmatch[[i]] <- 
        DF.Y[DF.Y$mz.Y > (mz.X[i] - 
                            (PPM/1e6*mz.X[i])) 
             & DF.Y$mz.Y < (mz.X[i] + 
                              (PPM/1e6*mz.X[i])) 
             & is.na(DF.Y$RT.Y), ]
      
    }
    Matched.X[i] <- as.numeric(nrow(MFmatch[[i]]))            
  }
  
  # Making a new data.frame to hold matched mass features
  DF.Y$Compound.Name.Y <- as.character(DF.Y$Compound.Name.Y)
  
  Matches <- data.frame(Compound.Name.Y = NA,
                        MassFeature.X = DF.X$MassFeature.X,
                        MassFeature.Y = NA, 
                        mz.X = DF.X$mz.X, 
                        mz.Y = NA,
                        RT.X = DF.X$RT.X,
                        RT.Y = NA,
                        NumMatched = Matched.X,
                        ppm = NA, 
                        RTdif = NA)
  
  Matches <- dlply(Matches, "MassFeature.X")
  
  # If there was more than 1 potential match and the RT was listed, take 
  # the one that was closest by m/z and then by RT. If there was no RT 
  # listed, take all matches. 
  # Calculate the difference in m/z in ppm and RT in min.
  for (i in names(Matched.X)[which(Matched.X > 0)]){
    
    for (n in 1:nrow(MFmatch[[i]])){
      MFmatch[[i]]$ppm[n] <- 
        abs((mz.X[i] - MFmatch[[i]]$mz.Y[n])/mz.X[i]*1e6)
      MFmatch[[i]]$RTdif[n] <- 
        abs(RT.X[i] - MFmatch[[i]]$RT.Y[n])
      MFmatch[[i]]$MassFeature.X <- i
    }
    MFmatch[[i]] <- arrange(MFmatch[[i]], ppm, RTdif)
    
    Matches[[i]]$MassFeature.Y <- 
      as.character(MFmatch[[i]]$MassFeature.Y[1])
    Matches[[i]]$mz.Y <- MFmatch[[i]]$mz.Y[1]
    Matches[[i]]$RT.Y <- MFmatch[[i]]$RT.Y[1]
    Matches[[i]]$ppm <- MFmatch[[i]]$ppm[1]
    Matches[[i]]$RTdif <- MFmatch[[i]]$RTdif[1]
    
  }
  
  #Matches is a simple df containing an index of Target Compounds and information on matching features in *.allpeaks.
  Matches <- rbind.fill(Matches)
  
  #Target.Compound.Areas is a df with a list of Target Compounds and corresponding peak areas
  names(Matches)[names(Matches) == 'MassFeature.Y'] <- 'MassFeature'
  
  Target.Compound.Areas <-
    Matches %>% join(xset.allpeaks, by = "MassFeature", type = "left")
  
  
  drops <-
    c(
      "MassFeature.X",
      "mz.X",
      "mz.Y",
      "mzmin",
      "mzmax",
      "RT.X",
      "rtmin",
      "rtmax",
      "NumMatched",
      "RT",
      "rt",
      "EddyQE_NoBins_XCMS"
    )
  Target.Compound.Areas <-
    Target.Compound.Areas[, !(names(Target.Compound.Areas) %in% drops)]
  
  names(Target.Compound.Areas)[names(Target.Compound.Areas) == 'RT.Y'] <-
    'RT'
  names(Target.Compound.Areas)[names(Target.Compound.Areas) == 'ppm'] <-
    'ppmdif'
  
  TargetCompoundsAreas_condensed <- Target.Compound.Areas[complete.cases(Target.Compound.Areas$mz),]
  
  #return(Matches)
  
# TCList <-
    #list(matches = Matches, tcl=Target.Compounds, tca = Target.Compound.Areas)
 #return(invisible(TCList))
setwd(ResultsDIR) 
write.csv(Matches, paste(Fraction, "Matches.csv", sep="."), row.names = FALSE) 
write.csv(Target.Compounds, paste(Fraction, "TargetCompounds.csv", sep="."), row.names = FALSE)  
write.csv(Target.Compound.Areas, paste(Fraction, "TargetCompoundAreas.csv", sep="."), row.names=FALSE)  

write.csv(TargetCompoundsAreas_condensed, paste(Fraction, "TargetCompoundAreasCondensed.csv", sep="."), row.names=FALSE)
 
  }
}

#Example:
#X <- Target.Compounds
#Y <- xsAq.allpeaks
#setwd("~/Desktop/Google Drive/3. Ingalls Lab/Compound info")
#source("Target.Compounds_finder.R")
#TCList <- mfmatch(X,Y,5,0.2)
#TCList$matches 
#TCList$tca
