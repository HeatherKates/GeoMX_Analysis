results <- c()
assayDataElement(object = target_Data.Q3norm, elt = "log_q") <-
  assayDataApply(target_Data.Q3norm, 2, FUN = log, base = 2, elt = "q_norm")
pData(target_Data.Q3norm)$testRegion <-factor(pData(target_Data.Q3norm)$pancreas_region)
target_Data.Q3norm@phenoData@varMetadata<- data.frame(colnames(target_Data.Q3norm@phenoData))

#make a subset of high 10% CV genes
# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_Data.Q3norm,
                         elt = "log_q", MARGIN = 1, calc_CV)
# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.9)]
subset <- subset(target_Data.Q3norm, TargetName %in% GOI)

pData(subset)$testRegion <-factor(pData(subset)$pancreas_region)
subset@phenoData@varMetadata<- data.frame(colnames(subset@phenoData))

#Create a dummy variable for mixed model
pData(subset)$Dummy <- rep(c("one","two","three","four","five","six","seven"),each=11)
pData(subset)$Dummy <-as.factor(pData(subset)$Dummy)

subset@phenoData@varMetadata<- data.frame(colnames(subset@phenoData))


for(ROI in c("islet_absent_ROI", "islet_present_ROI")) {
  
  ind <- pData(subset)$ROI_type == ROI
  mixedOutmc_pancregion <-
    mixedModelDE(subset[,ind],
                 elt = "log_q",
                 modelFormula = ~ testRegion+ (1 + testRegion | Dummy),
                 groupVar = "testRegion",
                 nCores = 1,
                 multiCore = FALSE,
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
  r_test$Subset <- ROI
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}
saveRDS(results,"P1.BetweenRegions.Result.RDATA")
