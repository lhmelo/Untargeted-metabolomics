# mfmatch function --------------------------------------------------

# This function will take two data.frames, NML and Y, and output a list containing two data frames. [[1]]  contains a listing which mass features match which in the other dataset. [[2]] contains peak areas accross all samples for the mass features listed in [[1]].

# Data.frame NML must include the following columns:
#       1. Compound.Name
#       2. MassFeature
#       3. mz
#       4. RT
#       
# List of target compounds should be data.frame NML:
# As of 27-01-17 using "MRM_Methods_Table.csv" for Cyano and HILIC targetted compounds and labelled standards
# 
# NMLsAq.allpeaks object from NMLCMS output should be data.frame Y.
# 
# Data.frame NML must include the following columns:
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



mfmatch <- function(PPM = 5, RTRange = 0.2){

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
    NML <- NeutralMassList #X
    
    
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
 
 
  TC <- Target.Compounds #Y
  
  DF.TC <- TC[, c("Compound.Name", "Compound.Type", "MassFeature", "mz", "RT")]
  DF.NML <- NML[, c("MassFeature", "mz", "RT")]
  
  DF.TC$MassFeature <- as.character(DF.TC$MassFeature)
  DF.NML$MassFeature <- as.character(DF.NML$MassFeature)
  
  DF.NML <- arrange(DF.NML, MassFeature)
  DF.TC <- arrange(DF.TC, MassFeature)
  
  # Renaming to keep track of which mass feature, m/z, and RT came from 
  # which dataset.
  names(DF.NML) <- paste(names(DF.NML), "NML", sep = ".")
  names(DF.TC) <- paste(names(DF.TC), "TC", sep = ".")
  
  MFmatch <- list()
  Matched.NML <- c()
  MFname.NML <- as.character(DF.NML$MassFeature.NML)
  mz.NML <- DF.NML$mz.NML
  names(mz.NML) <- DF.NML$MassFeature.NML
  RT.NML <- DF.NML$RT.NML
  names(RT.NML) <- DF.NML$MassFeature.NML
  
  # Checking each row in DF.NML for any matches in DF.Y.  The if statement is dealing with mass features for which RT is unknown.  This is unlikely to apply to my data.
  for (i in DF.NML$MassFeature.NML){
    
    print(i)
    MFmatch[[i]] <- DF.TC[
      DF.TC$mz.TC > (mz.NML[i] - (PPM/1e6*mz.NML[i])) 
      & DF.TC$mz.TC < (mz.NML[i] + (PPM/1e6*mz.NML[i])) 
      & DF.TC$RT.TC > RT.NML[i] - RTRange 
      & DF.TC$RT.TC < RT.NML[i] + RTRange, ]
    
   # if (nrow(MFmatch[[i]]) > 0) {
   #   MFmatch[[i]] <- 
   #     rbind(MFmatch[[i]], 
    #          DF.Y[DF.Y$mz.Y > (mz.NML[i] - 
   #                               (PPM/1e6*mz.NML[i])) 
   #                & DF.Y$mz.Y < (mz.NML[i] + 
   #                                 (PPM/1e6*mz.NML[i])) 
   #                & is.na(DF.Y$RT.Y), ])
      
      
   # } else {
    #  MFmatch[[i]] <- 
    #    DF.Y[DF.Y$mz.Y > (mz.NML[i] - 
   #                         (PPM/1e6*mz.NML[i])) 
    #         & DF.Y$mz.Y < (mz.NML[i] + 
    #                          (PPM/1e6*mz.NML[i])) 
    #         & is.na(DF.Y$RT.Y), ]
      
    #}
    Matched.NML[i] <- as.numeric(nrow(MFmatch[[i]]))            
  }
  
  # Making a new data.frame to hold matched mass features
  DF.NML$Compound.Name.NML <- as.character(DF.NML$Compound.Name.NML)
  
  Matches <- data.frame(Compound.Name.NML = DF.NML$Compound.Name.NML,
                        Compound.Type = DF.NML$Compound.Type.NML,
                        MassFeature.NML = DF.NML$MassFeature.NML,
                        MassFeature.Y = NA, 
                        mz.NML = DF.NML$mz.NML, 
                        mz.Y = NA,
                        RT.NML = DF.NML$RT.NML,
                        RT.Y = NA,
                        NumMatched = Matched.NML,
                        ppm = NA, 
                        RTdif = NA)
  
  #just splitting the data.frame, not applying a function.  Will rbind and recreate df after for loop.
  Matches <- dlply(Matches, "MassFeature.NML")
  
  # If there was more than 1 potential match and the RT was listed, take 
  # the one that was closest by m/z and then by RT. If there was no RT 
  # listed, take all matches. 
  # Calculate the difference in m/z in ppm and RT in min.
  for (i in names(Matched.NML)[which(Matched.NML > 0)]){
    
    for (n in 1:nrow(MFmatch[[i]])){
      MFmatch[[i]]$ppm[n] <- 
        abs((mz.NML[i] - MFmatch[[i]]$mz.Y[n])/mz.NML[i]*1e6)
      MFmatch[[i]]$RTdif[n] <- 
        abs(RT.NML[i] - MFmatch[[i]]$RT.Y[n])
      MFmatch[[i]]$MassFeature.NML <- i
    }
    MFmatch[[i]] <- arrange(MFmatch[[i]], ppm, RTdif)
    
    Matches[[i]]$MassFeature.Y <- 
      as.character(MFmatch[[i]]$MassFeature.Y[1])
    Matches[[i]]$mz.Y <- MFmatch[[i]]$mz.Y[1]
    Matches[[i]]$RT.Y <- MFmatch[[i]]$RT.Y[1]
    Matches[[i]]$ppm <- MFmatch[[i]]$ppm[1]
    Matches[[i]]$RTdif <- MFmatch[[i]]$RTdif[1]
    
  }
  
  #Matches is a simple df containing an indeNML of Target Compounds and information on matching features in *.allpeaks.
  Matches <- rbind.fill(Matches)
  
  #Target.Compound.Areas is a df with a list of Target Compounds and corresponding peak areas
  names(Matches)[names(Matches) == 'MassFeature.Y'] <- 'MassFeature'
  
  Target.Compound.Areas <-
    Matches %>% join(NMLset.allpeaks, by = "MassFeature", type = "left")
  
  
  drops <-
    c(
      "MassFeature.NML",
      "mz.NML",
      "mz.Y",
      "mzmin",
      "mzmaNML",
      "RT.NML",
      "rtmin",
      "rtmaNML",
      "NumMatched",
      "RT",
      "rt",
      "EddyQE_NoBins_NMLCMS"
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

#ENMLample:
#NML <- Target.Compounds
#Y <- NMLsAq.allpeaks
#setwd("~/Desktop/Google Drive/3. Ingalls Lab/Compound info")
#source("Target.Compounds_finder.R")
#TCList <- mfmatch(NML,Y,5,0.2)
#TCList$matches 
#TCList$tca
