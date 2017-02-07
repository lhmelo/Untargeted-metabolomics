# Untargeted-metabolomics

This repository contains code relavent to analysis of high resolution, full scan HPLC-MS data collected for untargeted metabolomic data. Examples of csv files required as inputs to run code (Params.csv and TargetCompoundList.csv) are in /data. 

Many functions have two versions.  An untagged version and one tagged with V2.  V2 functions are written to loop through multiple fractions of sample data (HILICNeg, HILICPos, CyanoAq and CyanoDCM) collected during a standard Ingalls Lab metabolomic workup. These functions can work on a single fraction as long as only one is specified in FractionList at the top of a (V2) workflow document.  The V2 versions are likely the most up to date versions of any given script.
