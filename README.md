# Untargeted-metabolomics

This repository contains code relavent to analysis of high resolution, full scan HPLC-MS data collected for untargeted metabolomic data. Examples of csv files required as inputs to run code (Params.csv and TargetCompoundList.csv) are in /data. 

Many functions have two versions.  An untagged version and one tagged with V2.  V2 functions are written to loop through multiple fractions of sample data (HILICNeg, HILICPos, CyanoAq and CyanoDCM) collected during a standard Ingalls Lab metabolomic workup. These functions can work on a single fraction as long as only one is specified in FractionList at the top of a (V2) workflow document.  The V2 versions are likely the most up to date versions of any given script.

Brief TOC:
(additional details in headers of individual scripts)

mfmaker_X.function.R

mfmatch.function.R

mfmatch_vesicles.R; mfmatch(); Takes as input two dataframes and outputs a dataframe containing only features which match in both dataframes.

mfNMLmatch_vesicles.function.R; mfNMLmatch(); Takes as input matcheddata from mfmatch_vesicles.function.R::mfmatch() and NML from NMLfunction() and outputs a dataframe detailing which mass features in matcheddaata might be related as adducts

NeutralMassList.function.R;  NMLfunction(); Takes as input xset.filtered and output from mfmaker. Produces a df of mass features that are signficant, by criteria assigned in mfmaker, and organized so that each row is a single mass feature which is related to other mass featuers as an adduct.

NMLmatch_test.V2.function.R: Dereplicates neutral mass list  produced by NMLfunction() using standard target compound list 
