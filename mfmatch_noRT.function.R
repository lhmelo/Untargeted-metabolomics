# mfmatch function --------------------------------------------------

# This function will take two data.frames, X and Y, and output a data.frame  
# listing which mass features match which in the other dataset. Data.frames 
# must include the following columns:
#       1. MassFeature
#       2. mz
#       3. RT

# DF.X should be xset.allpeaks
# DF.Y should be TargetCompoundWishList.csv located in "~/Desktop/MetabolomicsProjects/RScripts/data"
# Output is a data.frame of all mass features which match and what the m/z and RT differences are. The suffixes
# on the columns should make it clear which values came from which of the 
# original data.frames.
 
# If one of the datasets is much longer than the other, make the shorter one X. It will go faster.


mfmatch <- function(X, Y, PPM = 5){
  #require(plyr)
  require(stringr)
  
  DF.X <- X[, c("MassFeature", "mz", "RT")]
  DF.Y <- Y[, c("Compound.Name", "exact.mass")]
  
  DF.Y$MassFeature <- as.character(DF.Y$Compound.Name)
  DF.X$MassFeature <- as.character(DF.X$MassFeature)
  
  DF.X <- arrange(DF.X, MassFeature)
  DF.Y <- arrange(DF.Y, Compound.Name)
  
  # Renaming to keep track of which mass feature, m/z, and RT came from 
  # which dataset.
  

POLARITY <- ifelse(Params["POLARITY", Fraction]==1, "positive",
                    ifelse(Params["POLARITY", Fraction]==2, "positive",
                           ifelse(Params["POLARITY", Fraction]==3, "positive",
                                  ifelse(Params["POLARITY", Fraction]==4, "negative", NA))))
if (POLARITY == "positive") {
DF.Y$mz <- DF.Y$exact.mass + 1.0072
} else {
DF.Y$mz <- DF.Y$exact.mass - 1.0072
}

names(DF.X) <- paste(names(DF.X), "X", sep = ".")
names(DF.Y) <- paste(names(DF.Y), "Y", sep = ".")
  
  MFmatch <- list()
  Matched.X <- c()
  MFname.X <- as.character(DF.X$MassFeature)
  mz.X <- DF.X$mz
  names(mz.X) <- DF.X$MassFeature.X
  RT.X <- DF.X$RT
  names(RT.X) <- DF.X$MassFeature.X
  
  # Checking each row in DF.X for any matches in DF.Y
  for (i in DF.X$MassFeature.X){
    
    MFmatch[[i]] <- DF.Y[
      DF.Y$mz.Y > (mz.X[i] - (PPM/1e6*mz.X[i])) 
      & DF.Y$mz.Y < (mz.X[i] + (PPM/1e6*mz.X[i])),]
    
    Matched.X[i] <- as.numeric(nrow(MFmatch[[i]]))            
  }
  
  # Making a new data.frame to hold matched mass features
  Matches <- data.frame(MassFeature.X = DF.X$MassFeature.X,
                        Compound.Name.Y = NA, 
                        mz.X = DF.X$mz.X,
                        exact.mass.Y = NA,
                        mz.Y = NA,
                        RT.X = DF.X$RT.X,
                        RT.Y = NA,
                        NumMatched = Matched.X,
                        ppm = NA)
  
  Matches <- filter(Matches, NumMatched > 0)
  
  Matches <- dlply(Matches, "MassFeature.X")
  
  # If there was more than 1 potential match and the RT was listed, take 
  # the one that was closest by m/z and then by RT. If there was no RT 
  # listed, take all matches. 
  # Calculate the difference in m/z in ppm and RT in min.
  for (i in names(Matched.X)[which(Matched.X > 0)]){
    
    for (n in 1:nrow(MFmatch[[i]])){
      MFmatch[[i]]$ppm[n] <- 
        abs((mz.X[i] - MFmatch[[i]]$mz.Y[n])/mz.X[i]*1e6)
      MFmatch[[i]]$MassFeature.X <- i
    }
    MFmatch[[i]] <- arrange(MFmatch[[i]], ppm)
    
    Matches[[i]]$MassFeature.Y <- 
      as.character(MFmatch[[i]]$MassFeature.Y[1])
    Matches[[i]]$mz.Y <- MFmatch[[i]]$mz.Y[1]
    Matches[[i]]$exact.mass.Y <- MFmatch[[i]]$exact.mass.Y[1]
    Matches[[i]]$Compound.Name.Y <- MFmatch[[i]]$Compound.Name.Y[1]
    Matches[[i]]$RT.Y <- MFmatch[[i]]$RT.Y[1]
    Matches[[i]]$ppm <- MFmatch[[i]]$ppm[1]
    
    
  }
  
  Matches <- rbind.fill(Matches)
  
  Matches <- arrange(Matches, mz.X)
  
 # Matches <- Matches[complete.cases(Matches$mz.Y),]
  # 
  # names(Matches)
  # 
  # ColSuffix.X <- as.character(deparse(substitute(X)))
  # ColSuffix.Y <- as.character(deparse(substitute(Y)))
  # 
  # OrigNames <- names(Matches)
  # 
  # names(Matches)[str_detect(OrigNames, "Y$")] <- sub(
  #   "Y$", ColSuffix.Y, OrigNames)[str_detect(OrigNames, "Y$")]
  # names(Matches)[str_detect(OrigNames, "X$")] <- sub(
  #   "X$", ColSuffix.X, OrigNames)[str_detect(OrigNames, "X$")]

  return(Matches)
  
}
