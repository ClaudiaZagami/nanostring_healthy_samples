#install.packages("devtools")
#devtools::install_github("Nanostring-Biostats/NanoStringNCTools")
#devtools::install_github("Nanostring-Biostats/GeomxTools", ref = "dev")
#devtools::install_github("Nanostring-Biostats/GeoMxWorkflows", ref = "main")

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(readxl)

setwd("C:/Users/czagami/NanostringData")
datadir <- file.path("C:/Users/czagami/NanostringData/Data")
DCCdir <- file.path("C:/Users/czagami/NanostringData/Data/GeoMxNGS_Pipeline_Requeue_Normal_samples_CZ_DCC")

DCCFiles <- dir(file.path(datadir, "GeoMxNGS_Pipeline_Requeue_Normal_samples_CZ_DCC"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)

PKCFiles <- file.path(datadir, "Hs_R_NGS_WTA_v1.0.pkc")

SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

nano_healthy <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                       pkcFiles = PKCFiles,
                                       phenoDataFile = SampleAnnotationFile,
                                       phenoDataSheet = "Annotation template",
                                       phenoDataDccColName = "Sample_ID",
                                       protocolDataColNames = c("aoi", "roi"),
                                       experimentDataColNames = c("panel"))

#Asking which object is 
class(nano_healthy)

isS4(nano_healthy)
is(nano_healthy, "ExpressionSet")

# access the count matrix 
assayData(nano_healthy)[["exprs"]][1:3, 1:3] #genes and count per each DCC file
# access pheno data
pData(nano_healthy)[1:3, ]
# access the protocol data
pData(protocolData(nano_healthy))[1:3, ]
# access the probe information
fData(nano_healthy)[1:3, ]
# check feature type
featureType(nano_healthy)
# access PKC information
annotation(nano_healthy)

#in pheno data, there are three different annotations for E-cadherin
#make all E-cadherin the same annotation
#then, I want to add annotation for zone of the gland

library(stringr)
library(tidyverse)

#I created a vector with the right name and then substituted it to the original column.
#I was getting this error otherwise: Error in as.list.default(X) : 
#no method for coercing this S4 class to a vector

segment_names <- pData(nano_healthy)[, c("segment")]
segment_names <- str_replace_all(segment_names, "Cadherin", "E-Cadherin")
segment_names <- str_replace_all(segment_names, "E E-Cadherin", "E-Cadherin")
segment_names <- str_replace_all(segment_names, "E-E-Cadherin", "E-Cadherin")
segment_names


pData(nano_healthy)[, c("segment")] <- segment_names
pData(nano_healthy)[, c("segment")]

#Adding the pdata for the zones, using the same method

zone <- rep(c('Foveola','Foveola','Isthmus', 'Isthmus','Neck','Neck','Base','Base', 'Muscularis'), times=9)

pData(nano_healthy)[, c("Zone")] <- zone
pData(nano_healthy)[, c("Zone")]

#Adding an extra column for the simplified names of the patients. 

patient <- rep(c('27C', '29C', '31C'), each=27)
patient
pData(nano_healthy)[, c("patient")] <- patient

pData(nano_healthy)

#module used
library(knitr) #with this we can display tables
pkcs <- annotation(nano_healthy)
modules <- gsub(".pkc", "", pkcs) #function to replace strings
kable(data.frame(PKCs = pkcs, modules = modules))

#sample overview
library(dplyr) #data manipulation
library(ggforce) #extension of ggplot. 

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
#remember you called the column Zone with capital letter
count_mat <- count(pData(nano_healthy), `patient`, `segment`, `Zone`)

# gather the data and plot in order: class, slide name, region, segment
test_gr <- gather_set_data(count_mat, 1:3)
test_gr$x <- factor(test_gr$x,
                    levels = c("patient", "segment", "Zone"))

# plot Sankey
#this help to visualise the segmentation scheme and gives an overview of the samples distribution

ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = segment), alpha = 0.4, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(color = "white", size = 4) +
  theme_classic(base_size = 12) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") 

###QC and Pre-processing the data
# Shift counts to one
nano_healthy_shifted <- shiftCountsOne(nano_healthy, useDALogic = TRUE) #you do not want to use zero values

#Segment QC
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below

#Why did i change those values

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10) #this because we have low negative counts
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100) #in vimentin we have low nuclei. Nanostring specialst suggested to keep anything >20
       minArea = 1000)         # Minimum segment area (5000)


nano_healthy_shifted <-
  setSegmentQCFlags(nano_healthy_shifted, 
                    qcCutoffs = QC_params)   

# Collate QC Results
QCResults <- protocolData(nano_healthy_shifted)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))

QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

QC_Summary
#i have two flags in the aligned AOIs. 
#I will keep them for now anyway

#visualize segment QC
library(ggplot2)

col_by <- "segment"


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
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(nano_healthy_shifted), "Trimmed (%)", col_by, 80)
QC_histogram(sData(nano_healthy_shifted), "Stitched (%)", col_by, 80)
QC_histogram(sData(nano_healthy_shifted), "Aligned (%)", col_by, 80)
QC_histogram(sData(nano_healthy_shifted), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
#QC_histogram(nano_healthy_shifted, "nuclei", col_by, 20)

# calculate the negative geometric means for each module

negativeGeoMeans <- 
  esBy(negativeControlSubset(nano_healthy_shifted), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(nano_healthy_shifted)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData

negCols <- paste0("NegGeoMean_", modules) #this function crates a new part in my object. 
pData(nano_healthy_shifted)[, negCols] <- sData(nano_healthy_shifted)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(nano_healthy_shifted), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}
#there are different fragments that fall behind the dotted line. However, I am still ot sure about the meaning of the dotted line. 

# detatch neg_geomean columns ahead of aggregateCounts call
pData(nano_healthy_shifted) <- pData(nano_healthy_shifted)[, !colnames(pData(nano_healthy_shifted)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(nano_healthy_shifted)$NTC),
      col.names = c("NTC Count", "# of Segments"))

kable(QC_Summary, caption = "QC Summary Table for each Segment")


#do I need to remove the flagged segments or not? those are two vimentin segments 

dim(nano_healthy_shifted)


##PROBE QC
# Generally keep the qcCutoffs parameters unchanged. 
#Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers

nano_healthy_shift_QC <- setBioProbeQCFlags(nano_healthy_shifted, 
                             qcCutoffs = list(minProbeRatio = 0.1,
                                              percentFailGrubbs = 20), 
                             removeLocalOutliers = TRUE)

ProbeQCResults <- fData(nano_healthy_shift_QC)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

#exclude oulier probes

#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(nano_healthy_shift_QC, 
         fData(nano_healthy_shift_QC)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(nano_healthy_shift_QC)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

nano_healthy_shift_QC <- ProbeQCPassed 

##create gene-level count data

# Check how many unique targets the object has
length(unique(featureData(nano_healthy_shift_QC)[["TargetName"]]))

# collapse to targets
target_nano_healthy_SQC <- aggregateCounts(nano_healthy_shift_QC)
dim(target_nano_healthy_SQC)

exprs(target_nano_healthy_SQC)[1:5, 1:2]

###LOQ###
#limit of quantification LOQ based on the distribution of negative probes
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_nano_healthy_SQC))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_nano_healthy_SQC)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_nano_healthy_SQC)[, vars[1]] * 
             pData(target_nano_healthy_SQC)[, vars[2]] ^ cutoff)
  }
}
pData(target_nano_healthy_SQC)$LOQ <- LOQ

#Filtering out either segments and/or genes with abnormally low signal. 

#determining the number of genes detected in each segment across the dataset
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_nano_healthy_SQC)$Module == module
  Mat_i <- t(esApply(target_nano_healthy_SQC[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_nano_healthy_SQC)$TargetName, ]

#####segment gene detection 
#Filter out segments with exceptionally low signal. 
# Save detection rate information to pheno data

pData(target_nano_healthy_SQC)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_nano_healthy_SQC)$GeneDetectionRate <-
  pData(target_nano_healthy_SQC)$GenesDetected / nrow(target_nano_healthy_SQC)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_nano_healthy_SQC)$DetectionThreshold <- 
  cut(pData(target_nano_healthy_SQC)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_nano_healthy_SQC),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = pData(target_nano_healthy_SQC)$segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

ggplot(pData(target_nano_healthy_SQC),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = pData(target_nano_healthy_SQC)$patient)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "patients")

ggplot(pData(target_nano_healthy_SQC),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = pData(target_nano_healthy_SQC)$Zone)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Zones")

kable(table(pData(target_nano_healthy_SQC)$DetectionThreshold,
            pData(target_nano_healthy_SQC)$segment))

kable(table(pData(target_nano_healthy_SQC)$DetectionThreshold,
            pData(target_nano_healthy_SQC)$patient))

