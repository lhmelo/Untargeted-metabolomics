setwd("~/Desktop/EddyQE_XCMS")

library(IPO)

# Getting samples and peak picking -------------------------------------------------------
Samples <- list.files( getwd(), pattern=".mzXML", full.names=F, 
                       recursive = TRUE)   

Aq.Samples <- Samples[grepl(".Aq.", Samples)]
DCM.Samples <- Samples[grepl(".DCM.", Samples)]
Aq.Samples <- Aq.Samples[!grepl("Poo|Blk", Aq.Samples)]
DCM.Samples <- DCM.Samples[!grepl("Poo|Blk", DCM.Samples)]
Poo.Samples <- Samples[grepl("Poo_DCM", Samples)]

#Optimize peak picking parameters

peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$min_peakwidth <- c(4,10) 
peakpickingParameters$max_peakwidth <- c(10,20)
peakpickingParameters$ppm <- 5
resultPeakpicking <- 
  optimizeXcmsSet(files = Poo.Samples, 
                  params = peakpickingParameters, 
                  nSlaves = 2, 
                  subdir = 'rsmDirectory')
optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset
save(resultPeakpicking, file = "resultPeakpicking.RData")
save(optimizedXcmsSetObject, file="optimizedXcmsSetObject.RData")

#Optimize MORE peak picking parameters----

peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$min_peakwidth <- c(5,8) 
peakpickingParameters$max_peakwidth <- c(15,25)
peakpickingParameters$value_of_prefilter <- c(5000,10000)
peakpickingParameters$snthresh <- c(4,10)
peakpickingParameters$ppm <- 5
resultPeakpicking2 <- 
  optimizeXcmsSet(files = Poo.Samples, 
                  params = peakpickingParameters, 
                  nSlaves = 2, 
                  subdir = 'rsmDirectory')
optimizedXcmsSetObject2 <- resultPeakpicking2$best_settings$xset
save(resultPeakpicking2, file = "resultPeakpicking2.RData")
save(optimizedXcmsSetObject2, file="optimizedXcmsSetObject2.RData")

resultPeakpicking2$best_settings$parameters

#Use station 12 instead of Poo- since Poo were run as DDA at different resolution.

Station12.Samples <- list.files(path="~/Desktop/EddyQE_XCMS/Station12", pattern=".mzXML", full.names=T, recursive = TRUE)
Station12.Aq.Samples <- Station12.Samples[grepl("Aq", Station12.Samples)]
setwd("~/Desktop/EddyQE_XCMS")

peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$min_peakwidth <- c(5,7) 
peakpickingParameters$max_peakwidth <- c(15,25)
peakpickingParameters$value_of_prefilter <- c(10000,15000)
peakpickingParameters$snthresh <- c(4,10)
peakpickingParameters$ppm <- 5
peakpickingParameters$prefilter <- c(3,9)
resultPeakpicking3 <- 
  optimizeXcmsSet(files = Station12.Aq.Samples, 
                  params = peakpickingParameters, 
                  nSlaves = 2, 
                  subdir = 'rsmDirectory')
optimizedXcmsSetObject3 <- resultPeakpicking3$best_settings$xset
save(resultPeakpicking3, file = "resultPeakpicking3.RData")
save(optimizedXcmsSetObject3, file="optimizedXcmsSetObject3.RData")

resultPeakpicking3$best_settings$parameters

#Optimize retention time correction and grouping parameters----

retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
retcorGroupParameters$minfrac <- 0.5
retcorGroupParameters$mzwid <- c(0.001,0.01)
retcorGroupParameters$bw <- c(3,6)
retcorGroupParameters$max <- c(50,100)
resultRetcorGroup <-
  optimizeRetGroup(xset = optimizedXcmsSetObject3, 
                   params = retcorGroupParameters, 
                   nSlaves = 2, 
                   subdir = "rsmDirectory")



#Display optimized settings----

writeRScript(resultPeakpicking3$best_settings$parameters, 
             resultRetcorGroup$best_settings, 
             nSlaves=2)


xset <- xcmsSet(
  method = "centWave",
  peakwidth = c(5.9, 13.5), ppm = 5, noise = 0, snthresh =
    3.76, mzdiff = 0.000100000000000001, prefilter = c(1.9, 7000), mzCenterFun =
    "wMean", integrate = 1, fitgauss = FALSE, verbose.columns = FALSE, nSlaves =
    2
)

xset <- retcor(
  xset, method = "obiwarp",
  plottype = "none", distFunc = "cor_opt", profStep = 1, center =
    3, response = 1, gapInit = 0.78, gapExtend = 2.688, factorDiag = 2, factorGap =
    1, localAlignment = 0
)

xset <- group(
  xset, method = "density",
  bw = 0.25, mzwid = 0.00307, minfrac = 0.5, minsamp = 1, max =
    52.8
)

xset <- fillPeaks(xset, nSlaves = 2)