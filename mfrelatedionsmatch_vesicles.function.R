# mfNMLmatch function --------------------------------------------------

# This function will take two data.frames, X and Y, and output a data.frame  
# listing which mass features match which in the other dataset. Data.frames 
# must include the following columns:
#       1. MassFeature
#       2. mz
#       3. RT
#       4. log2ratio
#       5. AveTreat2
#       6. pvalue

# Output is a data.frame of all mass features which match and what the m/z and RT differences are. The suffixes
# on the columns should make it clear which values came from which of the 
# original data.frames.
 
# If one of the datasets is much longer than the other, make the shorter one X. It will go faster.


mfRImatch <- function(X, PPM = 5, RTRange = 0.2){
  require(plyr)
  require(stringr)
  
  ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
  setwd(ResultsDIR)
  load(paste(Fraction, "OtherIons", sep=""))
  # NeutralMassList <-
  #   read.csv(paste(Fraction, "NeutralMassList.csv", sep="."))
  
  
  Y <- OtherIons[!duplicated(OtherIons), ]
  
  Y_NM <- Y[complete.cases(Y$MassOfM),]
  
  Y_NM$UniqueMassFeature <- paste("I", round((Y_NM$mz),digits=4), 
                               "R", round( Y_NM$RT, digits=2), 
                               "M", round( Y_NM$MassOfM, digits=0), sep="")
  
  Y_Iso <- Y[complete.cases(Y$mzOfI),]
  
  Y_Iso$UniqueMassFeature <- paste("I", round((Y_Iso$mz),digits=4), 
                                 "R", round( Y_Iso$RT, digits=2), 
                                 "M", round( Y_Iso$mzOfI, digits=0), sep="")
  Y <- rbind(Y_NM, Y_Iso)
  
  
  DF.X <- X[, c("MassFeature", "log2ratio", "mz", "RT", "AveTreat2", "pvalue")]
  DF.Y <- Y[, c("MassFeature", "UniqueMassFeature", "mz", "RT", "IonType", "pcgroup", "MassOfM", "IsoGroup", "mzOfI")]
  
  DF.Y$MassFeature <- as.character(DF.Y$MassFeature)
  DF.X$MassFeature <- as.character(DF.X$MassFeature)
  DF.Y$IonType <- as.character(DF.Y$IonType)
  
  
  DF.X <- arrange(DF.X, MassFeature)
  DF.Y <- arrange(DF.Y, MassFeature)
  
  # Renaming to keep track of which mass feature, m/z, and RT came from 
  # which dataset.
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
      & DF.Y$mz.Y < (mz.X[i] + (PPM/1e6*mz.X[i])) 
      & DF.Y$RT.Y > RT.X[i] - RTRange 
      & DF.Y$RT.Y < RT.X[i] + RTRange, ]
    
    Matched.X[i] <- as.numeric(nrow(MFmatch[[i]]))            
  }
  
  
  
  # Making a new data.frame to hold matched mass features
  Matches <- data.frame(MassFeature.X = DF.X$MassFeature.X,
                        MassFeature.Y = NA, 
                        UniqueMassFeature.Y = NA,
                        IonType.Y = NA,
                        MassOfM.Y = NA,
                        mzofI.Y = NA,
                        mz.X = DF.X$mz.X, 
                        mz.Y = NA,
                        RT.X = DF.X$RT.X,
                        pcgroup.Y = NA,
                        IsoGroup.Y = NA,
                        AveTreat2.X = DF.X$AveTreat2.X,
                        log2ratio.X = DF.X$log2ratio.X,
                        pvalue.X = DF.X$pvalue.X,
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
    Matches[[i]]$UniqueMassFeature.Y <- 
      as.character(MFmatch[[i]]$UniqueMassFeature.Y[1])
    Matches[[i]]$mz.Y <- MFmatch[[i]]$mz.Y[1]
    Matches[[i]]$pcgroup.Y <- MFmatch[[i]]$pcgroup.Y[1]
    Matches[[i]]$IsoGroup.Y <- MFmatch[[i]]$IsoGroup.Y[1]
    Matches[[i]]$ppm <- MFmatch[[i]]$ppm[1]
    Matches[[i]]$RTdif <- MFmatch[[i]]$RTdif[1]
    Matches[[i]]$IonType.Y <- MFmatch[[i]]$IonType.Y[1]
    Matches[[i]]$MassOfM.Y <- MFmatch[[i]]$MassOfM.Y[1]
    Matches[[i]]$mzOfM.Y <- MFmatch[[i]]$mzOfM.Y[1]
    
  }
  
 
  
  
  Matches <- rbind.fill(Matches)
  
  # # Add to the Matches data.frame any MFs in Y that weren't found in X.
  # NoMatch <- DF.Y[DF.Y$MassFeature %in% 
  #                   setdiff(DF.Y$MassFeature, Matches$MassFeature.Y), ]
  
  #Matches <- rbind.fill(Matches, NoMatch)
  #Matches$mz.mean <- apply(Matches[, c("mz.X", "mz.Y")], 1, FUN = mean, 
  #                         na.rm=T)
  Matches <- arrange(Matches, RT.X)
  
 #Matches <- Matches[complete.cases(Matches$mz.Y),]
  
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
