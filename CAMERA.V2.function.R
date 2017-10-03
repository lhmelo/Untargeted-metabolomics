
#title: "Camera_Function"
#author: "Laura Hmelo"
#date: "December 19, 2016"

#This function processes an xcmsSet object using CAMERA to identify mass features which are related to one another as adducts or isotopomers.
#
#Input:
#1. a processed xset object from xcmsfunction.  Should be called xset3.
#
#Output:
#1. A list (CAMERAlist) containing [[1]] a processed peak list containing information about adducts, isotopes and neutral masses (xset.annot) and [[2]], a CAMERA object (xsaFA).

#Mode is "positive" or "negative". Default is positive


camera <- function(PPM = 5) {

library(CAMERA)
library(snow)
  
  for (j in 1:length(FractionList)){
    Fraction <- FractionList[j]
    print(Fraction)
    ResultsDIR <- as.character(Dirs[Fraction, "ResultsDIR"])
    setwd(ResultsDIR)
   load(paste(Fraction, "xset3.RData", sep="."))
   
   POLARITY <- ifelse(Params["POLARITY", Fraction]==1, "positive",
                      ifelse(Params["POLARITY", Fraction]==2, "positive",
                             ifelse(Params["POLARITY", Fraction]==3, "positive",
                                    ifelse(Params["POLARITY", Fraction]==4, "negative", NA))))
  
xsa <- xsAnnotate(xset3, nSlaves = 4)
print("xsAnnotate object created")

xsaF <- groupFWHM(xsa, perfwhm = 0.6)

# peak group verification with peakshape correlation (groupCorr); slow step----

xsaC <- groupCorr(xsaF)

#annotation of possible isotopes (findIsotopes)----
xsaFI <- findIsotopes(xsaC, ppm = PPM)

#annotation of adducts and calculating hypothetical masses for the group; very slow step---- 
#use xsaFA object to visualize EICs

xsaFA <- findAdducts(xsaFI, ppm = PPM, polarity = POLARITY)

# Make a table of the annotated peaks-------
xset.annot <- getPeaklist(xsaFA)

save(xsaFA, file=(paste(Fraction, "xsaFA", sep=".")))
save(xset.annot, file=(paste(Fraction, "xset.annot", sep=".")))

#CAMERAlist <- list(xset.annot, xsaFA)

#return(invisible(CAMERAlist))
print(paste("done", Fraction, sep=" "))

}
}