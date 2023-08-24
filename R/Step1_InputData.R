library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(knitr)
library(dplyr)
library(ggforce)
library(data.table)
library(stringr)
library(DT)

#Create the GeoMX object
Data <-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                      pkcFiles = PKCFiles,
                                      phenoDataFile = SampleAnnotationFile,
                                      phenoDataSheet = "Annotations",
                                      phenoDataDccColName = "Sample_ID",
                                      protocolDataColNames = c("aoi","roi"),
                                      experimentDataColNames = c("panel"))

pkcs <- annotation(Data)
modules <- gsub(".pkc", "", pkcs)
