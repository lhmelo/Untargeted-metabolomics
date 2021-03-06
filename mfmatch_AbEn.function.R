# mfmatch_AbEn function --------------------------------------------------

# This function will take two data.frames, X and Y, and output a data.frame  
# listing which mass features match which in the other dataset. Data.frames 
#  X should contain the Most Abundant Ions and RIs in vesicles, which is output from mfRImatch() (mfrelatedionsmatch_vesicles.function.R).  Data.frame Y is one of the logplot tables which are output from logplot.function.R.  Data.frames must contain columns listed in lines .. and .. below.  

# Output is a data.frame of all mass features which match and what the m/z and RT differences are. The suffixes
# on the columns should make it clear which values came from which of the 
# original data.frames.
 
# .X columns originate from xsets processed with both strains' vesicle data, whereas .Y columns originate from xsets processed with cell pellet and vesicle data.  MassFeatures and RTs may be slightly different due to slight RT corrections between the two groups.  
# 
# log2ratio and pvalue in final table refer to comparisons between cell pellets and vesicles.  Vesicles are Treat2.



mfmatch_AbEn <- function(X, Y, PPM = 5, RTRange = 0.2){
  require(plyr)
  require(stringr)
  
  DF.X <- X[, c("MassFeature", "mz", "RT", "IonType", "MassOfM", "mzOfI")]
  DF.Y <- Y[, c("MassFeature", "mz", "RT", "AveTreat1", "AveTreat2", "log2ratio", "pvalue")]
  
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
    
    Matched.X[i] <- as.numeric(nrow(MFmatch[[i]]))            
  }
  
  # Making a new data.frame to hold matched mass features
  Matches <- data.frame(MassFeature.X = DF.X$MassFeature.X,
                        MassFeature.Y = NA, 
                        mz.X = DF.X$mz.X, 
                        mz.Y = NA,
                        RT.X = DF.X$RT.X,
                        RT.Y = NA,
                        NumMatched = Matched.X,
                        IonType.X = DF.X$IonType.X,
                        MassOfM.X = DF.X$MassOfM.X,
                        mzOfI.X =DF.X$mzOfI.X,
                        AveTreat1.Y = NA,
                        AveTreat2.Y = NA,
                        log2ratio.Y = NA,
                        pvalue.Y = NA,
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
    Matches[[i]]$AveTreat1.Y <- MFmatch[[i]]$AveTreat1.Y[1]
    Matches[[i]]$AveTreat2.Y <- MFmatch[[i]]$AveTreat2.Y[1]
    Matches[[i]]$log2ratio.Y <- MFmatch[[i]]$log2ratio.Y[1]
    Matches[[i]]$pvalue.Y <- MFmatch[[i]]$pvalue.Y[1]
  }
  
  Matches <- rbind.fill(Matches)
  
  # # Add to the Matches data.frame any MFs in Y that weren't found in X.
  # NoMatch <- DF.Y[DF.Y$MassFeature %in% 
  #                   setdiff(DF.Y$MassFeature, Matches$MassFeature.Y), ]
  
  #Matches <- rbind.fill(Matches, NoMatch)
  #Matches$mz.mean <- apply(Matches[, c("mz.X", "mz.Y")], 1, FUN = mean, 
  #                         na.rm=T)
  Matches <- arrange(Matches, mz.X)
  
 Matches <- Matches[complete.cases(Matches$mz.Y),]
  
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
