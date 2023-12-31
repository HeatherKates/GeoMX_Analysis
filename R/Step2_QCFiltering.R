# Shift counts one
Data <- shiftCountsOne(Data, useDALogic = TRUE)

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

Data <-setSegmentQCFlags(Data, qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(Data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Samples, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}


all_samples <- rownames(pData(Data))

#Removal
Data.red <- Data[, QCResults$QCStatus == "PASS"]

Data.fail <- Data[, QCResults$QCStatus != "PASS"]

if(nrow(Data.fail@phenoData@data)>0){
failed.breakdown <- data.frame(table(Data.fail@phenoData@data$`slide name`,Data.fail@phenoData@data$segment,Data.fail@phenoData@data$ROI_type,Data.fail@phenoData@data$pancreas_region)) %>% filter(Freq>0)
colnames(failed.breakdown) <- c("Slide ID","Segmentation","ROI_type","Description","n segments")
datatable(failed.breakdown,caption = 'Breakdown of samples removed following QC)')
}

#Number of samples post-QC filtering
dim(Data.red)

#Probe-QC


# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
Data.red <- setBioProbeQCFlags(Data.red, 
                                   qcCutoffs = list(minProbeRatio = 0.1,
                                                    percentFailGrubbs = 20), 
                                   removeLocalOutliers = TRUE)

ProbeQCResults <- fData(Data.red)[["QCFlags"]]


#Make an object that is a subset of probes that passed QC
ProbeQCPassed <- 
  subset(Data.red, 
         fData(Data.red)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(Data.red)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)


#LOQ QC


# Check how many unique targets the object has
length(unique(featureData(Data)[["TargetName"]]))

# collapse to targets
target_Data <- aggregateCounts(ProbeQCPassed)
#dim(target_Data)

#
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_Data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_Data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_Data)[, vars[1]] * 
             pData(target_Data)[, vars[2]] ^ cutoff)
  }
}
pData(target_Data)$LOQ <- LOQ

#Filtering
#After determining the limit of quantification (LOQ) per segment, we recommend filtering out either segments and/or genes with #abnormally low signal. 
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_Data)$Module == module
  Mat_i <- t(esApply(target_Data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_Data)$TargetName, ]
#saveRDS(target_Data,"target_Data.RData")

#Segment gene detection
# Save detection rate information to pheno data
#target_Data <- readRDS("target_Data.RData")
pData(target_Data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_Data)$GeneDetectionRate <-
  pData(target_Data)$GenesDetected / nrow(target_Data)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_Data)$DetectionThreshold <- 
  cut(pData(target_Data)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

target_Data.filt <-
  target_Data[, pData(target_Data)$GeneDetectionRate >= .05]


library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_Data.filt)]
fData(target_Data.filt)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_Data.filt)$DetectionRate <-
  fData(target_Data.filt)$DetectedSegments / nrow(pData(target_Data.filt))

# Gene of interest detection table
##MARKER GENE HEATMAPS
#Markers retrieved from Azimuth reference
#https://azimuth.hubmapconsortium.org/references/human_pancreas/

CellType_markers <- read.csv("/blue/timgarrett/hkates/Campbell-Thompson/CellType_Markers.csv",header=FALSE)
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
  Number = fData(target_Data.filt)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_Data.filt)[goi, "DetectionRate"]),
  Label=marker_table %>% filter(gene==goi) %>% select(Expanded_Label)
)%>%distinct()


# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_Data.filt)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_Data.filt))
rownames(plot_detect) <- plot_detect$Freq

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_Data.filt), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_Data.filt2 <- 
  target_Data.filt[fData(target_Data.filt)$DetectionRate >= 0.1 |
                        fData(target_Data.filt)$TargetName %in% neg_probes, ]

# retain only detected genes of interest
goi.2 <- goi[goi %in% rownames(target_Data.filt2)]


library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "ROI_type"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_Data.filt2)),
             Segment = colnames(exprs(target_Data.filt2)),
             Annotation = pData(target_Data.filt2)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_Data.filt2), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_Data.filt2)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")


# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_Data.Q3norm <- GeomxTools::normalize(target_Data.filt2,
                                 norm_method = "quant", 
                                 desiredQuantile = .75,
                                 toElt = "q_norm")

