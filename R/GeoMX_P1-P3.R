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
library(readxl)
library(standR)
library(ggalluvial)
source("Share/geomx_import_data.HRK.R") #this replaces the standR sourcecode for these functions
#NOTES TO DO:
#IMAGES FROM INSTRUMENT, DOWNLOAD SCREENVIEW FOR EACH ROI
#DO SINGLE-CELL DECONVOLUTION FOR P1, P2, P3 using Azimuth/HPAP reference
#USE AZIMUTH GENE LISTS TO SHOW DE AMONG ISLET-PRESENT VS ISLET-ABSENT (marker_gene_table.csv)
#SHOW PCA ETC with POINTS ANNOTATED BY REGION ETC. (figure out if batch-correction is needed prior to)
#The only thing we can do with P1 is regional differences because only islets were captured from endocrine

#Read in all metadata data

P1 <- read_xlsx("Share/P1/P1_LabWorksheet.metadata.xlsx")
P2 <- read_xlsx("Share/P2/P2_LabWorksheet.metadata.xlsx")
P3 <- read_xlsx("Share/P3/P3_LabWorksheet.metadata.xlsx")

#Combine the metadata
Metadata <- rbind(P1,P2,P3)
#Drop the rows that were contrils
Metadata <- Metadata %>% dplyr::filter(!is.na(slide_name))
rownames(Metadata) <- Metadata$Sample_ID
Metadata <- Metadata %>% dplyr::rename("SegmentDisplayName"="Sample_ID")
Metadata$ROICoordinateX=0
Metadata$ROICoordinateY=0

#Read in all count data
P1_count <- read.csv("Share/P1/P1.countMatrix.csv")
P2_count <- read.csv("Share/P2/P2.countMatrix.csv")
P3_count <- read.csv("Share/P3/P3.countMatrix.csv")

#Combine all the count data
Count <- cbind(P1_count,P2_count,P3_count)
colnames(Count) <- gsub("\\.","-",colnames(Count)); colnames(Count) <- gsub("-dcc","",colnames(Count))
rownames(Count) <- Count$Gene
Count <- Count[,colnames(Count)%in%Metadata$SegmentDisplayName]
Count$TargetName <- rownames(Count)
Count <- Count %>% dplyr::relocate(TargetName)

#Feature Annotation
featureAnnoFile <- read.csv("featureFile.csv") %>% dplyr::rename("RTS_ID"="X")
rownames(featureAnnoFile) <- featureAnnoFile$RTS_ID

#Write files required by readGeoMX
saveRDS(Count,"P1-P3.count.RDATA")
saveRDS(Metadata,"P1-P3.metadata.RDATA")
saveRDS(featureAnnoFile,"featureAnnoFile.RDATA")



spe <- readGeoMx(Count, Metadata, featureAnnoFile,NegProbeName="NegProbe-WTX")# hasNegProbe = TRUE)

plotSampleInfo(spe, column2plot = c("ROI_type","pancreas_region","aoi_target"),textsize = 3)
spe <- addPerROIQC(spe, rm_genes = TRUE)
plotGeneQC(spe, ordannots = "aoi_target", col = "aoi_target", point_size = 2)
plotROIQC(spe, y_threshold = 50000, col = slide_name,
          x_axis = "area",
          y_axis = "lib_size",
          x_lab = "area",
          y_lab = "Library size")
spe <- spe[,rownames(colData(spe))[colData(spe)$lib_size > 50000]]
plotRLExpr(spe, ordannots = "slide_name", assay = 1, col = slide_name)
plotRLExpr(spe, ordannots = "slide_name", assay = 2, col = slide_name)

drawPCA(spe, assay = 2, col = slide_name, shape = aoi_target)
standR::plotMDS(spe, assay = 2, col = slide_name, shape = aoi_target_type)

saveRDS(spe,"spe_P1P2P3.RDATA")

#Variables
#datadir <- "/home/hkates/blue_garrett/Campbell-Thompson/GeoMX/Share/P2"
#DCCFiles <- dir(datadir, pattern=".*dcc$", full.names=TRUE)
#PKCFiles <- "Hs_R_NGS_WTA_v1.0.pkc"
#SampleAnnotationFile <- "P3_metadata.xlsx"
#phenoData <- "Metadata"
#marker_table <- read.csv("Other_Acinar_Endothelial_markers.csv",header=TRUE)

datadir <- "/home/hkates/blue_garrett/Campbell-Thompson/GeoMX/Share"
DCCFiles <- dir(datadir, pattern=".*1001660006002.*dcc$", full.names=TRUE)
PKCFiles <- "P1/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- "P1/P1_LabWorksheet.metadata.xlsx"
phenoData <- "Sheet1"
#marker_table <- read.csv("Other_Acinar_Endothelial_markers.csv",header=TRUE)

P1Data <-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                pkcFiles = PKCFiles,
                                phenoDataFile = SampleAnnotationFile,
                                phenoDataSheet = phenoData,
                                phenoDataDccColName = "Sample_ID",
                                #protocolDataColNames = c("aoi","roi"),
                                experimentDataColNames = c("panel"))


#To just write an output count matrix
# Shift counts one
GeoMxSet <- shiftCountsOne(P1Data, useDALogic = TRUE)
GeoMxSet <- aggregateCounts(GeoMxSet)
countData <- as.data.frame(GeoMxSet@assayData[["exprs"]])
countData$Gene <- rownames(GeoMxSet@assayData[["exprs"]])
countData <- countData %>% dplyr::relocate(Gene)
write.csv(file="P1/P1.countMatrix.csv",countData,row.names = FALSE)


#GeoMxSet <- subset(GeoMxSet,select=!is.na(phenoData(GeoMxSet)[["ROI_type"]])) #this is prob needed for an upstream error; check later


#write files for standR
write.csv(file="countData.csv",GeoMxSet@assayData[["exprs"]],row.names = TRUE)
write.csv(file="annoData.csv",GeoMxSet@phenoData@data,row.names = TRUE)
write.csv(file="featureFile.csv",GeoMxSet@featureData@data,row.names = TRUE)

##Select Segment QC
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)

GeoMxSet <-setSegmentQCFlags(GeoMxSet, qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(GeoMxSet)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

#Remove samples that failed QC
GeoMxSet <- GeoMxSet[, QCResults$QCStatus == "PASS"]
GeoMxSet <- subset(GeoMxSet,select=!is.na(phenoData(GeoMxSet)[["ROI_type"]])) #this is prob needed for an upstream error; check later

#Number of samples post-QC filtering
dim(GeoMxSet)

#Probe-QC

# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
GeoMxSet <- setBioProbeQCFlags(GeoMxSet, 
                             qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(GeoMxSet)[["QCFlags"]]


#Make an object that is a subset of probes that passed QC
GeoMxSet <- 
  subset(GeoMxSet, 
         fData(GeoMxSet)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(GeoMxSet)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)


#LOQ QC


# Check how many unique targets the object has
#length(unique(featureData(GeoMxSet)[["TargetName"]]))

# collapse to targets
GeoMxSet <- aggregateCounts(GeoMxSet)
#dim(target_GeoMxSet)

#
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
pkcs <- annotation(GeoMxSet)
modules <- gsub(".pkc", "", pkcs) 

LOQ <- data.frame(row.names = colnames(GeoMxSet))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(GeoMxSet)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(GeoMxSet)[, vars[1]] * 
             pData(GeoMxSet)[, vars[2]] ^ cutoff)
  }
}
pData(GeoMxSet)$LOQ <- LOQ

#Filtering
#After determining the limit of quantification (LOQ) per segment, we recommend filtering out either segments and/or genes with #abnormally low signal. 
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(GeoMxSet)$Module == module
  Mat_i <- t(esApply(GeoMxSet[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(GeoMxSet)$TargetName, ]
#saveRDS(GeoMxSet,"GeoMxSet.RData")

#Segment gene detection
# Save detection rate information to pheno data
#GeoMxSet <- readRDS("GeoMxSet.RData")
pData(GeoMxSet)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(GeoMxSet)$GeneDetectionRate <-
  pData(GeoMxSet)$GenesDetected / nrow(GeoMxSet)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(GeoMxSet)$DetectionThreshold <- 
  cut(pData(GeoMxSet)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))



GeoMxSet <-
  GeoMxSet[, pData(GeoMxSet)$GeneDetectionRate >= .05]


library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(GeoMxSet)]
fData(GeoMxSet)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(GeoMxSet)$DetectionRate <-
  fData(GeoMxSet)$DetectedSegments / nrow(pData(GeoMxSet))

# Gene of interest detection table
##MARKER GENE HEATMAPS
#Markers retrieved from Azimuth reference
#https://azimuth.hubmapconsortium.org/references/human_pancreas/

CellType_markers <- read.csv("../CellType_Markers.csv",header=FALSE)
colnames(CellType_markers) <- c("Label","Expanded_Label","OBO ontology ID","Markers")
CellType_markers$Markers <-  gsub(pattern = ";",replacement = ",",CellType_markers$Markers)
CellType_markers$Markers <-  gsub(pattern = " ",replacement = "",CellType_markers$Markers)
all_markers <- unlist(str_split(unlist(CellType_markers$Markers),pattern=","))
marker_table <- cbind(all_markers[1],CellType_markers[grepl(all_markers[1],CellType_markers$Markers),]) %>% dplyr::rename(gene = 1)
for (i in 2:length(all_markers)){
  marker_table <<- rbind(cbind(all_markers[i],CellType_markers[grepl(all_markers[i],CellType_markers$Markers),]) %>% dplyr::rename(gene = 1),marker_table)
}

goi <- marker_table$gene
goi_df <- data.frame(
  Gene = goi,
  Number = fData(GeoMxSet)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(GeoMxSet)[goi, "DetectionRate"]),
  Label=marker_table %>% filter(gene==goi) %>% dplyr::select(Expanded_Label)
)%>%distinct()


# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(GeoMxSet)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(GeoMxSet))
rownames(plot_detect) <- plot_detect$Freq

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(GeoMxSet), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
GeoMxSet <- 
  GeoMxSet[fData(GeoMxSet)$DetectionRate >= 0.1 |
                    fData(GeoMxSet)$TargetName %in% neg_probes, ]

# retain only detected genes of interest
goi <- goi[goi %in% rownames(GeoMxSet)]


library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "ROI_type"
Stat_data <- 
  data.frame(row.names = colnames(exprs(GeoMxSet)),
             Segment = colnames(exprs(GeoMxSet)),
             Annotation = pData(GeoMxSet)[, ann_of_interest],
             Q3 = unlist(apply(exprs(GeoMxSet), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(GeoMxSet)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")


# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
GeoMxSet <- GeomxTools::normalize(GeoMxSet,
                                                norm_method = "quant", 
                                                desiredQuantile = .75,
                                                fromElt="exprs",
                                                toElt = "q_norm")
#########
##STATS##
#########
results <- c()
assayDataElement(object = GeoMxSet , elt = "log_q") <-
  assayDataApply(GeoMxSet , 2, FUN = log, base = 2, elt = "q_norm")
pData(GeoMxSet)$testRegion <-factor(pData(GeoMxSet )$ROI_type)
pData(GeoMxSet)$aoi_Target <-factor(pData(GeoMxSet )$aoi_target)
pData(GeoMxSet)$SlideID <-factor(pData(GeoMxSet )$slide_name)

#Subset to omit immune type ROI
subset <- subset(GeoMxSet ,select=phenoData(GeoMxSet )[["ROI_type"]] %in% c("islet_present_ROI","islet_absent_ROI"))

#marker_table %>% filter(!aoi_target==marker) %>% dplyr::select(gene) %>% unlist()
#For each segment type, create a subset of the Geomxset object by first removing marker genes for all other segment types
#(use marker_table)
marker_table <- read.csv("Other_Acinar_Endothelial_markers.csv",header=TRUE)
marker_table$Genes <- paste(marker_table$Gene,marker_table$Gene_synonym)
marker_table$Genes <- gsub(";",",",marker_table$Genes)
marker_table$Genes <- gsub(" ",",",marker_table$Genes)
marker_table$Genes <- gsub(",,",",",marker_table$Genes)


for(marker in c("CD31+CD34+:endothelial_cells","None:Other","PanCK+:duct_cells")) {
  genes_marker <- gsub("+","",marker)
  #subset segments of the aoi type
  ind <- pData(subset)$aoi_Target == marker
  ##remove genes that are markers for other aoi types
  #make a list of marker genes for aoi types other than the aoi of interest
  keep_genes <- marker_table %>% filter(grepl(genes_marker,Tissue)) %>% dplyr::select(Genes) %>% unlist()
  #code to subset remove the marker genes from the data -used in function below
  #subset(subset[, ind], !TargetName %in% drop_genes)
  
  mixedOutmc_pancregion <-
    mixedModelDE(subset(subset[, ind], TargetName  %in% unlist(str_split(string = keep_genes,pattern = ","))),
                 elt = "log_q",
                 modelFormula = ~ testRegion + (1 | SlideID),
                 groupVar = "testRegion",
                 nCores = 4,
                 multiCore = TRUE,
                 pAdjust = "BY",
                 pairwise=TRUE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc_pancregion["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc_pancregion),
                  rep, nrow(mixedOutmc_pancregion["lsmeans", ][[1]])))
  r_test$Subset <- marker
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}

colnames(results) <- c("Gene", "Subset", "contrast", "log2FoldChange", 
                       "pvalue", "padj")
data="P3_3292023"
saveRDS(results,paste0(data,".EachAOItype_betweenROItype.Result.RDATA"))
saveRDS(subset,paste0(data,".NanoStringGeoMXSet.postQC.NoImmune.RDATA"))


