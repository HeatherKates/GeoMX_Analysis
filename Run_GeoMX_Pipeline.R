#Read in data files
Case="P1"
datadir <- "Example_Data/GeoMx_NGS_Pipeline_02_04_2022_5_05_33-523906740/GeoMx_NGS_Pipeline_02_04_2022_5_05_33-ds.9555984bb7004944a9f4a6ff1ff8c256/"
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- "Example_Data/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- "Example_Data/GeoMX_P1_Annotations.xlsx"

source("R/Step1_InputData.R")
source("R/Step2_QCFiltering.R")
source("R/Step3_LMM.R")
saveRDS(target_Data.Q3norm,file="P1.QCfiltered.RDATA")