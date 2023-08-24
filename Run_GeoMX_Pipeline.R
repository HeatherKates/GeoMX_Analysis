#Read in data files
datadir <- "Example_Data/GeoMx_NGS_Pipeline_02_04_2022_5_05_33-523906740/"
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- "Example_Data/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- "Example_Data/"