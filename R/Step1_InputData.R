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

#Read in data files
datadir <- "/orange/timgarrett/GE5701-MThompson-S1-HW525DRXY-Lane1-327318995/GeoMx_NGS_Pipeline_02_04_2022_5_05_33-523906740/GeoMx_NGS_Pipeline_02_04_2022_5_05_33-ds.9555984bb7004944a9f4a6ff1ff8c256/"
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- "/blue/timgarrett/hkates/Campbell-Thompson/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- "/blue/timgarrett/hkates/Campbell-Thompson/PancDB/CThompsonAnnotations.xlsx"

#Create the GeoMX object
ThompsonData <-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                      pkcFiles = PKCFiles,
                                      phenoDataFile = SampleAnnotationFile,
                                      phenoDataSheet = "Annotations",
                                      phenoDataDccColName = "Sample_ID",
                                      protocolDataColNames = c("aoi","roi"),
                                      experimentDataColNames = c("panel"))

pkcs <- annotation(ThompsonData)
modules <- gsub(".pkc", "", pkcs)

#simplify slide name
ThompsonData@phenoData@data$`slide name` <- gsub("WTA 12-15-21","",ThompsonData@phenoData@data$`slide name`)
ThompsonData@phenoData@data <- ThompsonData@phenoData@data %>% mutate(ROI_type=recode(compartment, endocrine = "islet_present_ROI", exocrine = "islet_absent_ROI"))
#ThompsonData@phenoData@data <- ThompsonData@phenoData@data %>% mutate(Marker_Target=recode(region, acinar = "None:Other" , vessel = "CD31+CD34+:endothelial_cells", beta_cells="INS+:beta_cells",duct="PanCK+:duct_cells",immune_other="None:Other"))
#ThompsonData@phenoData@data <- ThompsonData@phenoData@data %>% mutate(Target=recode(region, acinar = "Other" , vessel = "endothelial_cells", beta_cells="beta_cells",duct="duct_cells",immune_other="Other"))
#Pancreas region
ThompsonData@phenoData@data <- ThompsonData@phenoData@data %>% rename(pancreas_region=organ_region)
#Make these data objects match
ThompsonData@phenoData@varMetadata<- data.frame(colnames(ThompsonData@phenoData))
