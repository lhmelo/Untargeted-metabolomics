#title: "NeutralMassList.function"
#author: "Laura Hmelo"
#last update: "January 20, 2017"

# This script takes the data.frame containing filtered xcms output (xset.filtered) and 
# 1) filters out only ions which are significant by the metric chosen in mfmaker (e.g., ions with distributions p<0.05 by a statistical test) and 
# 2) further refines this list by of significant ions by seperating out ions with true neutral masses as well as calculating 'hypothetical' neutral masses for ions which do not associate with other ions as adducts.  
# 
# 
# NMLfunction----
#
#Inputs:
#1: sig.groupnames: a data frame containing columns "groupname" and "Pvalue"; this is output from the XCMS script
#2: *.filtered: a data frame containing columns "groupname", "MassFeature", "AvgSmp"; this dataframe is output from XCMS script.
#3: *otherions: a data.frame containing columns MassFeature, MZ, RT, pcgroup, IonType,Charge, MassofM, IsoGroup, MassofI.  This data.frame is output from function camera in CAMERA.function.R
#
#Outputs is a list, NMLists, which contains:
#[[1]] UniqueNML: A one column dataframe containg unique neutral masses
#[[2]] UniqueIsoHypNM: A one column dataframe containing unique values of hypothetical neutral masses calculated from ions which were associated with another isotope, but not an adduct.  This assumes the ions are M+H
#[[3]] UniqueHypNM: A one column dataframe containing unique values of hypothetical neutral masses calculated from ions which did not associate with any other ions as adducts or istopes.  This assumes all ions are M+H.
#[[4]] NeutralMassList: A dataframe which contains all ions which associate with another ion as an adduct.  Each neutral mass is listed as a seperate line but a single line should not be considered a single metabolite. A single metabolite may be represented by several adducts. This dataframe retains information about peak area in each sample.

NMLfunction <- function(OtherIons, xset.filtered, sig.groupnames) {

require(plyr)
require (dplyr)

#df sig.groupnames will define the feature list moving forward.  For now, it represents features which have passed QC and have pvalue < 0.5 by ANOVA.  This excludes any features represented evenly accross all samples.
#This line pulls any adduct and isotope identifications found by CAMERA and matching one of the significant mass features
#Many features have no Adduct or Isotope information.  No adducts or isotopes were associated with these features. I will deal with these later


OtherIons <- OtherIons[, !(colnames(OtherIons) %in% c("mz"))]

#SigFeatures <- sig.groupnames %>% join(xset.filtered, by="groupname", type="left") %>% 
#  select(MassFeature, groupname, mz, AveSmp, Pvalue) %>%
#  join(OtherIons, by="MassFeature", type="left")
#  
SigFeatures <- sig.groupnames %>% join(xset.filtered, by="groupname", type="left") %>% 
  select(-Fvalue, -mzmin, -mzmax, -rt, -rtmin, -rtmax, -npeaks, -RT, -AveSmp) %>%
  join(OtherIons, by="MassFeature", type="left")

#generate a list of Significant features which associated with an adduct or isotope ID

SigFeatures_wAddIso <- SigFeatures[complete.cases(SigFeatures$pcgroup),]

#Only list features for which CAMERA was able to identify a NM

NeutralMassList <- SigFeatures[complete.cases(SigFeatures$MassOfM),]

#Remove redundant listings and generate a csv file which can be used to input into Metlin.

UniqueNML <- as.data.frame(unique(NeutralMassList$MassOfM))


#generate a list of significant features which have an isotope, but are not already binned in with adduct list
#
IsoList <- SigFeatures[complete.cases(SigFeatures$mzOfI),]
IsoList$HypNM <- IsoList$mzOfI - 1.0078
IsoList <- IsoList[which(IsoList[,'IonType'] == "I"),]
UniqueIsoHypNM <- as.data.frame(unique(IsoList$HypNM))

#generate a list of Significant features which DID NOT associate with an adduct or isotope ID

SigFeatures_woAddIso <- SigFeatures[!complete.cases(SigFeatures$pcgroup),]

#Make a new column for a hypothesized Neutral Mass, subtract a hydrogen

SigFeatures_woAddIso$HypNM <- SigFeatures_woAddIso$mz  - 1.0078
UniqueHypNM <-  as.data.frame(unique(SigFeatures_woAddIso$HypNM))

NMLists <- list(UniqueNML, UniqueIsoHypNM, UniqueHypNM, NeutralMassList)
return(invisible(NMLists))

}
