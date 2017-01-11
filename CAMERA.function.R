
#title: "Camera_Function"
#author: "Laura Hmelo"
#date: "December 19, 2016"

#Mode is "positive" or "negative". Default is positive

setwd(outputpath)

camera <- function(xset, Polarity ="positive", PPM = 5) {

library(CAMERA)
library(snow)

xsa <- xsAnnotate(xset)

xsaF <- groupFWHM(xsa, perfwhm = 0.6)

# peak group verification with peakshape correlation (groupCorr); slow step----

xsaC <- groupCorr(xsaF)

#annotation of possible isotopes (findIsotopes)----
xsaFI <- findIsotopes(xsaC, ppm = PPM)

#annotation of adducts and calculating hypothetical masses for the group; very slow step---- 
#use xsaFA object to visualize EICs

xsaFA <- findAdducts(xsaFI, ppm = PPM, polarity = Polarity)

# Make a table of the annotated peaks-------
xset.annot <- getPeaklist(xsaFA)

CAMERAlist <- list(xset.annot, xsaFA)

return(invisible(CAMERAlist))

}
