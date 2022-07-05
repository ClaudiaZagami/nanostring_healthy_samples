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

# Asking which object is 
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

# in pheno data, there are three different annotations for E-cadherin
# make all E-cadherin the same annotation
# then, I want to add annotation for zone of the gland

library(stringr)
library(tidyverse)

# I created a vector with the right name and then substituted it to the original column.
# I was getting this error otherwise: Error in as.list.default(X) : 
# no method for coercing this S4 class to a vector

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

# Adding an extra column for epithelium vs stroma only. 
# in this way we can do comparisons among the two classes only if we want to go easier 
# we can also find clusters only upon those characteristics

class <- rep(c('epithelium', 'stroma', 'stroma', 'epithelium', 'stroma', 'epithelium', 'stroma', 'epithelium', 'stroma' ), times=9)
class
pData(nano_healthy)[, c("class")] <- class

pData(nano_healthy)

# Adding extra column for nuclei count

nuclei <- c(393, 114, 217, 362,448, 765, 571, 718, 267, 435, 85, 151, 510, 218, 787, 178, 617, 400, 348,75, 216, 317, 265, 599, 229, 738, 347, 406, 113, 55, 227, 155, 442, 128, 382, 221, 346, 43, 41, 183, 133, 597, 134, 391, 176, 351, 79, 68, 367, 116, 703, 110, 500, 355, 525, 26, 181, 416, 36, 115, 75, 451, 256, 395, 60, 63, 213, 95, 567, 79, 439, 130, 16, 2, 41, 89, 46, 118, 98, 558, 201)
pData(nano_healthy)[, c("nuclei")] <- nuclei

pData(nano_healthy)

# module used
library(knitr) # with this we can display tables
pkcs <- annotation(nano_healthy)
modules <- gsub(".pkc", "", pkcs) # function to replace strings
kable(data.frame(PKCs = pkcs, modules = modules))

# sample overview
library(dplyr) # data manipulation
library(ggforce) # extension of ggplot. 

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
# remember you called the column Zone with capital letter
count_mat <- count(pData(nano_healthy), `patient`, `class`, `segment`, `Zone`, `nuclei`)

# gather the data and plot in order: class, slide name, region, segment
test_gr <- gather_set_data(count_mat, 1:4)
test_gr$x <- factor(test_gr$x,
                    levels = c("patient", "class", "segment", "Zone"))

# plot Sankey
# this help to visualise the segmentation scheme and gives an overview of the samples distribution

ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = segment), alpha = 0.4, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(color = "white", size = 4) +
  theme_classic(base_size = 12) + 
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") 

### QC and Pre-processing the data
# Shift counts to one
nano_healthy_shifted <- shiftCountsOne(nano_healthy, useDALogic = TRUE) #you do not want to use zero values

# Segment QC
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below

# Why did i change those values

QC_params <-
  list(
    minSegmentReads = 1000, # Minimum number of reads (1000)
    percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
    percentStitched = 80,   # Minimum % of reads stitched (80%)
    percentAligned = 80,    # Minimum % of reads aligned (80%)
    percentSaturation = 50, # Minimum sequencing saturation (50%)
    minNegativeCount = 1,   # Minimum negative control counts (10) #this because we have low negative counts
    maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
    minNuclei = 20,         # Minimum # of nuclei estimated (100) #in vimentin we have low nuclei. Nanostring specialst suggested to keep anything >20
    minArea = 1000
    )         # Minimum segment area (5000)


nano_healthy_shifted <-
  setSegmentQCFlags(nano_healthy_shifted, 
    qcCutoffs = QC_params
    )   

# Collate QC Results
QCResults <- protocolData(nano_healthy_shifted)[["QCFlags"]]

flag_columns <- colnames(QCResults)

QC_Summary <- data.frame(
  Pass = colSums(!QCResults[, flag_columns]),
  Warning = colSums(QCResults[, flag_columns])
  )

QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

QC_Summary["TOTAL FLAGS", ] <-
  c(
    sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING")
    )

QC_Summary

# i have two flags in the aligned AOIs. 
# I will keep them for now anyway

# visualize segment QC
library(ggplot2)

col_by <- "class"


QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(
    assay_data,
    aes_string(
      x = paste0("unlist(`", annotation, "`)"),
      fill = fill_by
      )
    ) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept = thr, lty = "dashed", color = "black") +
      theme_bw() + 
      guides(fill = "none") +
      facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
      labs(x = annotation, y = "Segments, #", title = annotation)
    if (!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}


QC_histogram(sData(nano_healthy_shifted), "Trimmed (%)", col_by, 80)
QC_histogram(sData(nano_healthy_shifted), "Stitched (%)", col_by, 80)
QC_histogram(sData(nano_healthy_shifted), "Aligned (%)", col_by, 80)
QC_histogram(sData(nano_healthy_shifted), "Saturated (%)", col_by, 50) +
  labs(
    title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)"
    )
QC_histogram(nano_healthy_shifted, "nuclei", col_by, 20)

# calculate the negative geometric means for each module

negativeGeoMeans <- 
  esBy(negativeControlSubset(nano_healthy_shifted), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }
      ) 
protocolData(nano_healthy_shifted)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData

negCols <- paste0("NegGeoMean_", modules) # this function crates a new part in my object. 
pData(nano_healthy_shifted)[, negCols] <- sData(nano_healthy_shifted)[["NegGeoMean"]]
for (ann in negCols) {
  plt <- QC_histogram(pData(nano_healthy_shifted), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

# there are different fragments that fall behind the dotted line. However, I am still not sure about the meaning of the dotted line. 

# detatch neg_geomean columns ahead of aggregateCounts call
pData(nano_healthy_shifted) <- pData(nano_healthy_shifted)[, !colnames(pData(nano_healthy_shifted)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(nano_healthy_shifted)$NTC),
  col.names = c("NTC Count", "# of Segments")
  )

kable(QC_Summary, caption = "QC Summary Table for each Segment")


# do I need to remove the flagged segments or not? those are two vimentin segments 

dim(nano_healthy_shifted)


## PROBE QC
# Generally keep the qcCutoffs parameters unchanged. 
# Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers

nano_healthy_shift_QC <- setBioProbeQCFlags(nano_healthy_shifted, 
  qcCutoffs = list(
    minProbeRatio = 0.1,
    percentFailGrubbs = 20
    ), 
    removeLocalOutliers = TRUE
  )

ProbeQCResults <- fData(nano_healthy_shift_QC)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(
  Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
  Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
  Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0 & 
    !ProbeQCResults$GlobalGrubbsOutlier)
)

#exclude outlier probes

# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(
    nano_healthy_shift_QC, 
    fData(nano_healthy_shift_QC)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
      fData(nano_healthy_shift_QC)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE
  )

dim(ProbeQCPassed)

nano_healthy_shift_QC <- ProbeQCPassed 

## create gene-level count data

# Check how many unique targets the object has
length(unique(featureData(nano_healthy_shift_QC)[["TargetName"]]))

# collapse to targets
target_nano_healthy_SQC <- aggregateCounts(nano_healthy_shift_QC)
dim(target_nano_healthy_SQC)

exprs(target_nano_healthy_SQC)[1:5, 1:2]

### LOQ###
# limit of quantification LOQ based on the distribution of negative probes
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_nano_healthy_SQC))
for(module in modules) {
  vars <- paste0(
    c("NegGeoMean_", "NegGeoSD_"),
    module
  )
  if(all(vars[1:2] %in% colnames(pData(target_nano_healthy_SQC)))) {
    LOQ[, module] <-
      pmax(
        minLOQ,
        pData(target_nano_healthy_SQC)[, vars[1]] * 
          pData(target_nano_healthy_SQC)[, vars[2]] ^ cutoff
      )
  }
}
pData(target_nano_healthy_SQC)$LOQ <- LOQ

# Filtering out either segments and/or genes with abnormally low signal. 

# determining the number of genes detected in each segment across the dataset
LOQ_Mat <- c()
for (module in modules) {
  ind <- fData(target_nano_healthy_SQC)$Module == module
  Mat_i <- t(esApply(target_nano_healthy_SQC[ind, ], 
    MARGIN = 1,
    FUN = function(x) {
      x > LOQ[, module]
    }
    ))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_nano_healthy_SQC)$TargetName, ]

##### segment gene detection 
# Filter out segments with exceptionally low signal. 
# Save detection rate information to pheno data

pData(target_nano_healthy_SQC)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_nano_healthy_SQC)$GeneDetectionRate <-
  pData(target_nano_healthy_SQC)$GenesDetected / nrow(target_nano_healthy_SQC)

###### How can I see for example 3%?
# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_nano_healthy_SQC)$DetectionThreshold <- 
  cut(pData(target_nano_healthy_SQC)$GeneDetectionRate,
    breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
    labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%")
  )

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(
  pData(target_nano_healthy_SQC),
  aes(x = DetectionThreshold)
) +
  geom_bar(aes(fill = pData(target_nano_healthy_SQC)$segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Gene Detection Rate",
    y = "Segments, #",
    fill = "Segment Type"
  )

ggplot(
  pData(target_nano_healthy_SQC),
  aes(x = DetectionThreshold)
) +
  geom_bar(aes(fill = pData(target_nano_healthy_SQC)$class)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Gene Detection Rate",
    y = "Segments, #",
    fill = "class"
  )

ggplot(
  pData(target_nano_healthy_SQC),
  aes(x = DetectionThreshold)
) +
  geom_bar(aes(fill = pData(target_nano_healthy_SQC)$patient)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Gene Detection Rate",
    y = "Segments, #",
    fill = "patients"
  )

ggplot(
  pData(target_nano_healthy_SQC),
  aes(x = DetectionThreshold)
) +
  geom_bar(aes(fill = pData(target_nano_healthy_SQC)$Zone)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Gene Detection Rate",
    y = "Segments, #",
    fill = "Zones"
  )

kable(table(
  pData(target_nano_healthy_SQC)$DetectionThreshold,
  pData(target_nano_healthy_SQC)$segment
))

kable(table(
  pData(target_nano_healthy_SQC)$DetectionThreshold,
  pData(target_nano_healthy_SQC)$patient
))

kable(table(
  pData(target_nano_healthy_SQC)$DetectionThreshold,
  pData(target_nano_healthy_SQC)$class
))

# I cannot filter the segments with 1-5% detection rate because I would delete all the data

## now I can determine the gene detection rate for genes across the study

library(scales) # for percent

# Can't really understand this bit of the code
# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_nano_healthy_SQC)]
fData(target_nano_healthy_SQC)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_nano_healthy_SQC)$DetectionRate <-
  fData(target_nano_healthy_SQC)$DetectedSegments / nrow(pData(target_nano_healthy_SQC))

# Gene of interest detection table
# this is just to understand the concept

genesoi <- c(
  "WNT8A", "WNT10A", "BMP8A", "MAPK6", "MAPKAPK3", "MAPK14",
  "PGC", "COL1A2", "COL5A2", "MUC5AC", "MUC6"
)
genesoi_df <- data.frame(
  Gene = genesoi,
  Number = fData(target_nano_healthy_SQC)[genesoi, "DetectedSegments"],
  DetectionRate = percent(fData(target_nano_healthy_SQC)[genesoi, "DetectionRate"])
)

# not really understanding the df I get from this function. 

## Gene filtering 

plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(
    c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
    function(x) {
      sum(fData(target_nano_healthy_SQC)$DetectionRate >= x)
    }
  ))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_nano_healthy_SQC))
rownames(plot_detect) <- plot_detect$Freq

# visualize
# total number of genes detected in different percentages of segments
# understand global gene detection in our study and select how many low detected 
# genes to filter out of the dataset

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
    vjust = 1.6, color = "black", size = 4
  ) +
  scale_fill_gradient2(
    low = "orange2", mid = "lightblue",
    high = "dodgerblue3", midpoint = 0.65,
    limits = c(0,1),
    labels = scales::percent
  ) +
  theme_bw() +
  scale_y_continuous(
    labels = scales::percent, limits = c(0,1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "% of Segments",
    y = "Genes Detected, % of Panel > LOQ"
  )

# 5% of the segments have 2062 genes detected so the thrashold I would use to 
# normalize has to be at least 5% otherwise, with 10% I would loose too many genes

# Subset to target genes detected in at least 5% of the samples.
# Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_nano_healthy_SQC), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_nano_healthy_norm <- 
  target_nano_healthy_SQC[fData(target_nano_healthy_SQC)$DetectionRate >= 0.05 |
    fData(target_nano_healthy_SQC)$TargetName %in% neg_probes, ]
dim(target_nano_healthy_norm)

#### NORMALISATION 
# The two common methods for normalization of DSP-NGS RNA data are 
# i) quartile 3 (Q3) or ii) background normalization
# Given the low negative probe counts we will use Q3. 
# there should be a separation between these two values to ensure we have stable measure of Q3 signal.

library(reshape2)  # for melt
library(cowplot) # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives

exprs(target_nano_healthy_norm)

ann_of_interest <- "class"
stat_norm <- 
  data.frame(
    row.names = colnames(exprs(target_nano_healthy_norm)),
    Segment = colnames(exprs(target_nano_healthy_norm)),
    Annotation = pData(nano_healthy)[,ann_of_interest],
    Q3 = unlist(apply(exprs(target_nano_healthy_norm), 2,
      quantile, 0.75, 
      na.rm = TRUE
    )),
    NegProbe = exprs(target_nano_healthy_norm)[neg_probes, ]
  )

Stat_norm_m <- melt(stat_norm, 
  measure.vars = c("Q3", "NegProbe"),
  variable.name = "Statistic", value.name = "Value"
)

plt1 <- ggplot(
  Stat_norm_m,
  aes(x = Value, fill = Statistic)
) +
  geom_histogram(bins = 40) + 
  theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt1

plt2 <- ggplot(
  stat_norm,
  aes(x = NegProbe, y = Q3, color = Annotation)
) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + 
  guides(color = "none") + 
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt2 

plt3 <- ggplot(
  stat_norm,
  aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)
) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + 
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
plt3

btm_row <- plot_grid(plt2, plt3, 
  nrow = 1, labels = c("B", ""),
  rel_widths = c(0.43,0.57)
)
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

#### NORMALISATION####
# I can create two dataset, one normalized Q3 and another one for background

#### Q3 norm (75th percentile) for WTA  with or without custom spike-ins
target_nano_healthy_nQ3 <- NanoStringNCTools::normalize(target_nano_healthy_norm , 
  data_type = "RNA",
  norm_method = "quant", 
  desiredQuantile = .75,
  toElt = "q_norm"
)

# Background normalization for WTA/CTA without custom spike-in
target_nano_healthy_backg <- normalize(target_nano_healthy_norm,
  data_type = "RNA",
  norm_method = "neg",
  fromElt = "exprs",
  toElt = "neg_norm"
)

# visualize the first 15 segments with each normalization method
boxplot(exprs(target_nano_healthy_nQ3)[,1:15],
  col = "#9EDAE5", main = "Raw Counts",
  log = "y", names = 1:15, xlab = "Segment",
  ylab = "Counts, Raw"
)


boxplot(assayDataElement(target_nano_healthy_nQ3[,1:15], elt = "q_norm"),
 col = "#2CA02C", main = "Q3 Norm Counts",
 log = "y", names = 1:15, xlab = "Segment",
 ylab = "Counts, Q3 Normalized"
 )

boxplot(assayDataElement(target_nano_healthy_backg[,1:15], elt = "neg_norm"),
 col = "#FF7F0E", main = "Neg Norm Counts",
 log = "y", names = 1:15, xlab = "Segment",
 ylab = "Counts, Neg. Normalized"
 )

# the plots are slightly different but I am not sure which one is the best to use. 


### unsupervised analysis 
#install.packages('umap')
library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
# what does happen if I add or remove seeds?
# Marton: You will get a slightly different output every time. Selecting a seed will
# allow you to be consistent within a run at least.

custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

umap_out <-
  umap(t(log2(assayDataElement(target_nano_healthy_nQ3 , elt = "q_norm"))),  
    config = custom_umap
  )
pData(target_nano_healthy_nQ3)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

Q3_umap <- ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = UMAP1, y = UMAP2, color = Zone, shape = segment)
) +
  geom_point(size = 3) +
  theme_bw()
Q3_umap

# umap q3 normalisation zone and class

ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = UMAP1, y = UMAP2, color = class, shape = Zone)
) +
  geom_point(size = 3) +
  theme_bw()

#change seeds to 50 to see what changes

custom_umap1 <- umap::umap.defaults
custom_umap1$random_state <- 50

umap_out1 <-
  umap(t(log2(assayDataElement(target_nano_healthy_nQ3 , elt = "q_norm"))),  
    config = custom_umap1
  )
pData(target_nano_healthy_nQ3)[, c("UMAP1", "UMAP2")] <- umap_out1$layout[, c(1,2)]
ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = UMAP1, y = UMAP2, color = Zone, shape = segment)) +
  geom_point(size = 3) +
  theme_bw()  # graph changes completely.

# need to understand which are the best conditions. 

# try with 3 components and 42 seeds

custom_umap2 <- umap::umap.defaults
custom_umap2$random_state <- 42
custom_umap2$n_components <- 3

umap_out2 <-
  umap(t(log2(assayDataElement(target_nano_healthy_nQ3 , elt = "q_norm"))),  
    config = custom_umap2
  )
pData(target_nano_healthy_nQ3)[, c("UMAP1", "UMAP2", "UMAP3")] <- umap_out2$layout[, c(1,2,3)]
ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = UMAP1, y = UMAP2, z = UMAP3, color = Zone, shape = segment)) +
  geom_point(size = 3) +
  theme_bw()

# Do not know how to add the z axes to the plot #however with the same number of seeds I got a different type of graph, why?

# now with the background normalisation 

umap_out_BG <-
  umap(t(log2(assayDataElement(target_nano_healthy_backg , elt = "neg_norm"))),  
    config = custom_umap
  )
pData(target_nano_healthy_backg)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

Bkg_umap <- ggplot(
  pData(target_nano_healthy_backg),
  aes(x = UMAP1, y = UMAP2, color = Zone, shape = segment)
) +
  geom_point(size = 3) +
  theme_bw()

Bkg_umap

# plot the umap from the two normalisations one next to the other to see the difference

library(ggpubr)

ggarrange(Bkg_umap + rremove("legend"), Q3_umap, labels = c("background norm", "Q3 norm"), ncol = 2, nrow = 1)
#to improve the position of the lables, and the distribution of the graphs. 

#it seems that the two normalisation methods do not give particular differences in the segregation of the samples. 
# Marton: correct!

# run tSNE
#seed 42, on Q3 normalised data and background. 
#plotting segments and zones
# Marton: not sure tSNE is needed if you already did UMAP

set.seed(42) # set the seed for tSNE as well

#Q3 normalisation 
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_nano_healthy_nQ3 , elt = "q_norm"))),
   perplexity = ncol(target_nano_healthy_nQ3)*.15
  )
pData(target_nano_healthy_nQ3)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]

Q3_tsne <- ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = tSNE1, y = tSNE2, color = segment, shape = zone)
) +
  geom_point(size = 3) +
  theme_bw()


# same for the normalised data upon background noise

tsne_out1 <-
  Rtsne(t(log2(assayDataElement(target_nano_healthy_backg , elt = "neg_norm"))),
   perplexity = ncol(target_nano_healthy_backg)*.15
  )
pData(target_nano_healthy_backg)[, c("tSNE1", "tSNE2")] <- tsne_out1$Y[, c(1,2)]

Bkg_tsne <- ggplot(
  pData(target_nano_healthy_backg),
  aes(x = tSNE1, y = tSNE2, color = segment, shape = zone)
) +
  geom_point(size = 3) +
  theme_bw()

# putting all the graph together. 

ggarrange(Bkg_umap + rremove("legend"), Q3_umap, Bkg_tsne + rremove("legend") , Q3_tsne,
  labels = c("background norm", "Q3 norm"), ncol = 2, nrow = 2
)

# to be improved: removing the legend and putting it as an extra so that all the graphs will appear equal. 
# the two graphs have a differnet color annotation, so I cannot use the same legend for all the graphs. 

#### CLUSTERING HIGH CV GENES####

# from here we might not need it since we might use the counts for our network analysis

library(pheatmap)  # for heatmap
# create a log2 transform of the data for analysis
# I will use the Q3 data for now

assayDataElement(object = target_nano_healthy_nQ3, elt = "log_q") <-
 assayDataApply(target_nano_healthy_nQ3, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {
  sd(x) / mean(x)
}
CV_Q3norm_healthy <- assayDataApply(target_nano_healthy_nQ3,
  elt = "log_q", MARGIN = 1, calc_CV
)
# show the highest CD genes and their CV values

sort(CV_Q3norm_healthy, decreasing = TRUE)[1:50]

# Identify genes in the top 3rd of the CV values

top3rd_Q3 <- names(CV_Q3norm_healthy)[CV_Q3norm_healthy > quantile(CV_Q3norm_healthy, 0.8)]
pheatmap(assayDataElement(target_nano_healthy_nQ3[top3rd_Q3, ], elt = "log_q"),
  scale = "row", 
  show_rownames = FALSE, show_colnames = FALSE,
  border_color = NA,
  clustering_method = "complete", #I tryied to change the method but it seems that the clusters do not make much sense. 
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  breaks = seq(-3, 3, 0.05),
  color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
  annotation_col = 
  pData(target_nano_healthy_nQ3)[, c("Zone", "segment", "patient", "class")]
) 
#how can I change the colours of the annotation?
# Marton: I changed it to complete linkage.
# Complete linkage clustering means that the distance from one cluster to another 
# is calculated based on the furthest members of the cluster. 
# The used clustering is sensitive for the furthest elements. 
# Complete linkage does not join together with the furthest clusters, 
# producing a clear picture.

### DIFFERENTIAL EXPRESSION

# within slide analysis
# taking in consideration zones: study of differences between morphological structures
### comparing structures that co-exist within the a given tissue we will use the LMM model with a random slope####
# Morphological structure (zone) is out test variable. 
# control for tissue sub sampling using slide name or patient
# Benjamin-Hochberg multiple test correction.

# convert test variables to factors
pData(target_nano_healthy_nQ3)$testClass <- 
  factor(pData(target_nano_healthy_nQ3)$class, c("epithelium", "stroma"))

pData(target_nano_healthy_nQ3)[["patient"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["patient"]])

pData(target_nano_healthy_nQ3)[["class"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["class"]])

pData(target_nano_healthy_nQ3)[["segment"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["segment"]])

pData(target_nano_healthy_nQ3)[["Zone"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["Zone"]])

assayDataElement(object = target_nano_healthy_nQ3, elt = "log_q") <-
  assayDataApply(target_nano_healthy_nQ3, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package

#I am not sure I wrote the formula correctly. I need to understand the terms I need to insert in the formula. 

# =====
# Marton: To understand the loop, let's do it for one example at first

ind <- pData(target_nano_healthy_nQ3)$Zone == "Foveola"
mixedOutmc <-
  mixedModelDE(target_nano_healthy_nQ3[, ind],
               elt = "log_q",
               modelFormula = ~ testClass + (1 + testClass | patient),
               groupVar = "testClass",
               nCores = parallel::detectCores(),
               multiCore = F)

# formatting the result to be more readable
r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(r_test)
r_test <- as.data.frame(r_test)
r_test$Contrast <- tests
r_test$Gene <- 
  unlist(lapply(colnames(mixedOutmc),
                rep, nrow(mixedOutmc["lsmeans", ][[1]])))
r_test$Zone <- "Foveola" # note that I renamed this variable
r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
r_test <- r_test[, c("Gene", "Zone", "Contrast", "Estimate", 
                     "Pr(>|t|)", "FDR")]

# =====

# Marton: Now do the loop again
# Marton: We create a list (which is a dictionary in other languages) to save the dataframes into

# (Marton: I took out Muscularis for now as it throws an error since it's not found in the epithelial layer)
# Marton: Differential expression between epithelium and stroma 
# for each morphological structure

results_1 <- list()

for (i in c("Foveola", "Isthmus", "Neck", "Base")) {
  
  ind <- pData(target_nano_healthy_nQ3)$Zone == i
    
  mixedOutmc <-
    mixedModelDE(target_nano_healthy_nQ3[, ind],
                 elt = "log_q",
                 modelFormula = ~ testClass + (1 + testClass | patient),
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = F
    )

  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Zone <- paste0(i) # note that I renamed this variable
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Zone", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  
  results_1[[paste0(i)]] <- r_test
}

# Now we can browse results_1 like so: results_1[["Foveola"]] or results_1[["Base"]] etc.
# Or put them in a single large dataframe, like so: 
  
all_results_from_epi_stroma <- bind_rows(results_1)
all_results_from_epi_stroma_1 <- bind_rows(results_1)

write_excel_csv(results_1[["Foveola"]], file = "results_foveola_EvsS.csv", ",")
write_excel_csv(results_1[["Isthmus"]], file = "results_isthmus_EvsS.csv", ",")
write_excel_csv(results_1[["Neck"]], file = "results_neck_EvsS.csv", ",")
write_excel_csv(results_1[["Base"]], file = "results_base_EvsS.csv", ",")

# Keep significant only
# REMEMBER: play with threshold
all_results_from_epi_stroma_significant <- all_results_from_epi_stroma |> dplyr::filter(abs(Estimate) >= 1 & FDR <= 0.5)
# Too losen, it returns 140 genes

all_results_from_epi_stroma_significant1 <- all_results_from_epi_stroma |> dplyr::filter(abs(Estimate) >= 0.25 & FDR <= 0.001)
# Parameters suggested by Marton but return only 1 gene

all_results_from_epi_stroma_significant4 <- all_results_from_epi_stroma |> dplyr::filter(abs(Estimate) >= 0.6 & `Pr(>|t|)` < 0.05)

all_results_from_epi_stroma_significant2 <- all_results_from_epi_stroma |> dplyr::filter(abs(Estimate) >= 1 & FDR <= 0.05)
#23 genes

all_results_from_epi_stroma_significant3 <- all_results_from_epi_stroma |> dplyr::filter(abs(Estimate) >= 0.5 & FDR <= 0.05)

# Results are all upregulated - the downregulated genes are not significant
# Plot below show significantly upregulated genes

ggplot(all_results_from_epi_stroma_significant, 
       aes(x = Gene, y = Estimate, fill=Zone)) +
  geom_bar(stat='identity')+
  coord_flip()+
  facet_wrap(~Zone)

# convert test variables to factors
pData(target_nano_healthy_nQ3)$testClass1 <- 
  factor(pData(target_nano_healthy_nQ3)$class, c("stroma", "epithelium"))

pData(target_nano_healthy_nQ3)[["patient"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["patient"]])

pData(target_nano_healthy_nQ3)[["class"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["class"]])

pData(target_nano_healthy_nQ3)[["segment"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["segment"]])

pData(target_nano_healthy_nQ3)[["Zone"]] <- 
  factor(pData(target_nano_healthy_nQ3)[["Zone"]])

assayDataElement(object = target_nano_healthy_nQ3, elt = "log_q") <-
  assayDataApply(target_nano_healthy_nQ3, 2, FUN = log, base = 2, elt = "q_norm")

results_2 <- list()

for (i in c("Foveola", "Isthmus", "Neck", "Base")) {
  
  ind <- pData(target_nano_healthy_nQ3)$Zone == i
  
  mixedOutmc <-
    mixedModelDE(target_nano_healthy_nQ3[, ind],
                 elt = "log_q",
                 modelFormula = ~ testClass1 + (1 + testClass1 | patient),
                 groupVar = "testClass1",
                 nCores = parallel::detectCores(),
                 multiCore = F
    )
  
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Zone <- paste0(i) # note that I renamed this variable
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Zone", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  
  results_2[[paste0(i)]] <- r_test
}


# Categorize Results based on P-value & FDR for plotting

all_results_from_epi_stroma_1$Color <- "NS or FC < 0.5"
all_results_from_epi_stroma_1$Color[all_results_from_epi_stroma_1$`Pr(>|t|)` < 0.05] <- "P < 0.05"
all_results_from_epi_stroma_1$Color[all_results_from_epi_stroma_1$FDR < 0.05] <- "FDR < 0.05"
all_results_from_epi_stroma_1$Color[all_results_from_epi_stroma_1$FDR < 0.001] <- "FDR < 0.001"
all_results_from_epi_stroma_1$Color[abs(all_results_from_epi_stroma_1$Estimate) < 0.5] <- "NS or FC < 0.5"
all_results_from_epi_stroma_1$Color <- factor(all_results_from_epi_stroma_1$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

all_results_from_epi_stroma_filt1 <- all_results_from_epi_stroma |> dplyr::filter(abs(Color) == P < 0.05 & FDR < 0.05 & FDR < 0.001 )

# volcano plot
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
all_results_from_epi_stroma$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" in our case epithelium
all_results_from_epi_stroma$diffexpressed[all_results_from_epi_stroma$Estimate > 0.6 & all_results_from_epi_stroma$`Pr(>|t|)` < 0.05] <- "Epithelium"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN" in our case stroma
all_results_from_epi_stroma$diffexpressed[all_results_from_epi_stroma$Estimate < -0.6 & all_results_from_epi_stroma$`Pr(>|t|)` < 0.05] <- "Stroma"

#cut off estimate 0.6 and FDR 0.05
res_epistr_1_filt <- all_results_from_epi_stroma |> dplyr::filter(abs(Estimate) >= 0.6 & (Estimate) <= -0.6 & FDR <= 0.05)
# with this params we get 0 results. 
#res_epistr_1_filt <- res_epistr_1_filt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" in our case epithelium
#res_epistr_1_filt$diffexpressed[res_epistr_1_filt$Estimate > 0.6 & res_epistr_1_filt$`Pr(>|t|)` < 0.05] <- "Epithelium"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN" in our case stroma
#res_epistr_1_filt$diffexpressed[res_epistr_1_filt$Estimate < -0.6 & res_epistr_1_filt$`Pr(>|t|)` < 0.05] <- "Stroma"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=all_results_from_epi_stroma, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed)) + geom_point() + theme_minimal()
p

# Add lines 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

## Change point color 

# by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("red", "black", "blue"))
p3

# create a vector to define the colours with a rule
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Stroma", "Epithelium", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

# Now write down the name of genes beside the points.
# Create a new column "delabel" to my detaframe, that will contain the name of genes differentially expressed (NA in case they are not)
all_results_from_epi_stroma$delabel <- NA
all_results_from_epi_stroma$delabel[all_results_from_epi_stroma$diffexpressed != "NO"] <- all_results_from_epi_stroma$Gene[all_results_from_epi_stroma$diffexpressed != "NO"]

ggplot(data=all_results_from_epi_stroma, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=all_results_from_epi_stroma, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "log2 fold change", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# ============filtering and pathway analysis for FOLEVOLA ====================
#filter for foveola

epi_stroma_foveola <- all_results_from_epi_stroma |> dplyr::filter(Zone == "Foveola")
epi_stroma_foveola

# creating df of only differntially expressed genes in foveola
epi_stroma_foveola_sig <- epi_stroma_foveola |> dplyr::filter(delabel != "NA")

#vulcano plot foveola only 

ggplot(data=epi_stroma_foveola, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "log2 fold change", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# Filtering significant in stroma
# cut off?
epi_stroma_foveola_sig1 <- epi_stroma_foveola |> dplyr::filter(abs(Estimate) >= 0.6 & FDR <= 0.05)

#create excel file 
write_excel_csv(epi_stroma_foveola_sig, file = "epi_stroma_foveola_sig.csv", ",")


# create an excel file from results1
#write_excel_csv(all_results_from_epi_stroma, file = "results_1.csv", ",")

# Violin plot for example gene (WNT8A), using all data

ggplot(pData(target_nano_healthy_nQ3),
       aes(x = class, fill = class,
           y = assayDataElement(target_nano_healthy_nQ3["WNT8A", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "WNT8A Expression") +
  scale_y_continuous(trans = "log2") +
  facet_wrap(~Zone) + #check how to dispose it in order as the gland
  theme_bw()

#hes7
ggplot(pData(target_nano_healthy_nQ3),
       aes(x = class, fill = class,
           y = assayDataElement(target_nano_healthy_nQ3["HES7", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "HES7 Expression") +
  scale_y_continuous(trans = "log2") +
  facet_wrap(~Zone) +
  theme_bw()

# pathway analysis with ClusterProfiler for FOVEOLA ONLY

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)

colnames(epi_stroma_foveola_sig) <- c("Gene", "Zone", "Contrast", "log2FoldChange", "pvalue", "FDR", "diffexpressed", "delabel")
epi_stroma_foveola_sig

# SET THE DESIRED ORGANISM 
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

## Example to understand the code
# Taken from https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# reading in data from deseq2
#df = read.csv("drosphila_example_de.csv", header=TRUE)
# we want the log2 fold change 
#original_gene_list <- df$log2FoldChange
# name the vector
#names(original_gene_list) <- df$X
#original_gene_list
# omit any NA values 
#gene_list<-na.omit(original_gene_list)
#gene_list
# sort the list in decreasing order (required for clusterProfiler)
#gene_list = sort(gene_list, decreasing = TRUE)
#gene_list

# do the same with my data 
# log2 fold change
epi_stroma_fov_genelist <- epi_stroma_foveola_sig$log2FoldChange

# name the vector
names(epi_stroma_fov_genelist) <- epi_stroma_foveola_sig$Gene

# sort the list in decreasing order (required for clusterProfiler)
epi_stroma_fov_genelist = sort(epi_stroma_fov_genelist, decreasing = TRUE)
epi_stroma_fov_genelist

#check which options of keytype fo rthe gene names re available
keytypes(org.Hs.eg.db) # I will try using symbol

gse <- gseGO(epi_stroma_fov_genelist, 
                          ont ="ALL", 
                          keyType = "SYMBOL", 
                          nPerm = 10000, 
                          minGSSize = 3, 
                          maxGSSize = 800, 
                          pvalueCutoff = 0.05, 
                          verbose = TRUE, 
                          OrgDb = organism, 
                          pAdjustMethod = "none")

results_GSA_epistr_fov <- gse@result
write_excel_csv(results_GSA_epistr_fov, file = "results_GSA_epistr_fov.csv", ",")

# Dotplot
require(DOSE)
dotplot(gse, showCategory=30, split=".sign", label_format = 5) + facet_grid(.~.sign)

# Encrichment Map
emapplot(gse, showCategory = 10) 
# Error: Error in has_pairsim(x) : 
# Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.

#install.packages("ggnewscale")
library(ggnewscale)
ema <- pairwise_termsim(gse, method = "JC", semData = NULL, showCategory = 200)
emapplot(ema, showCategory = 10)

# Ridgeplot
#install.packages("ggridges")
library(ggridges)
ridgeplot(gse) + labs(x = "enrichment distribution")

cnetplot(gse, categorySize="pvalue", foldChange=epi_stroma_fov_genelist, showCategory = 3)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[3], geneSetID = 1)


#Insert bonferroni correction
gse_epistr_fov <- gseGO(epi_stroma_fov_genelist, 
                        ont ="ALL", 
                        keyType = "SYMBOL", 
                        nPerm = 10000, 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = organism, 
                        pAdjustMethod = "bonferroni")

#with bonferroni correction, here are no results. why?

# Create gseKEGG object
# Prepare input

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(epi_stroma_fov_genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
ids
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
dedup_ids
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
epi_stroma_foveola_sig2 = epi_stroma_foveola_sig[epi_stroma_foveola_sig$Gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
epi_stroma_foveola_sig2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_epistro_foveola <- epi_stroma_foveola_sig2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_epistro_foveola) <- epi_stroma_foveola_sig2$Y
kegg_epistro_foveola

# sort the list in decreasing order (required for clusterProfiler)
kegg_epistro_foveola = sort(kegg_epistro_foveola, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_epistro_foveola,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
# kegg is not working: error: Error in download.KEGG.Path(species) : 
# 'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...

# ============filtering and pathway analysis for isthmus ====================

#filter for isthmus

epi_stroma_isthmus <- all_results_from_epi_stroma |> dplyr::filter(Zone == "Isthmus")
epi_stroma_isthmus

# creating df of only differntially expressed genes in isthmus
epi_stroma_isthmus_sig <- epi_stroma_isthmus |> dplyr::filter(delabel != "NA")

#vulcano plot isthmus only 

ggplot(data=epi_stroma_isthmus, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "log2 fold change", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


#create excel file 
library(readr)
write_excel_csv(epi_stroma_isthmus_sig, file = "epi_stroma_isthmus_sig.csv", ",")

# pathway analysis with ClusterProfiler for ISTHMUS ONLY

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)

colnames(epi_stroma_isthmus_sig) <- c("Gene", "Zone", "Contrast", "log2FoldChange", "pvalue", "FDR", "diffexpressed", "delabel")
epi_stroma_isthmus_sig

# SET THE DESIRED ORGANISM 
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# log2 fold change
epi_stroma_ist_genelist <- epi_stroma_isthmus_sig$log2FoldChange

# name the vector
names(epi_stroma_ist_genelist) <- epi_stroma_isthmus_sig$Gene

# sort the list in decreasing order (required for clusterProfiler)
epi_stroma_ist_genelist = sort(epi_stroma_ist_genelist, decreasing = TRUE)
epi_stroma_ist_genelist

gse_epistr_ist <- gseGO(epi_stroma_ist_genelist, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

results_GSA_epistr_ist <- gse_epistr_ist@result
write_excel_csv(results_GSA_epistr_ist, file = "results_GSA_epistr_ist.csv", ",")
gse_epistr_ist@result

# Dotplot
require(DOSE)
dotplot(gse_epistr_ist, showCategory=30, split=".sign", label_format = 5) + facet_grid(.~.sign)

# Encrichment Map
emapplot(gse_epistr_ist, showCategory = 10) 
# Error: Error in has_pairsim(x) : 
# Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.

#install.packages("ggnewscale")
library(ggnewscale)
ema_ist <- pairwise_termsim(gse_epistr_ist, method = "JC", semData = NULL, showCategory = 200)
emapplot(ema_ist, showCategory = 10)

# Ridgeplot
#install.packages("ggridges")
library(ggridges)
ridgeplot(gse_epistr_ist) + labs(x = "enrichment distribution")

cnetplot(gse_epistr_ist, categorySize="pvalue", foldChange=epi_stroma_fov_genelist, showCategory = 3)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse_epistr_ist, by = "all", title = gse_epistr_ist$Description[3], geneSetID = 1)

# ============filtering and pathway analysis for NECK ====================

#filter for neck

epi_stroma_neck <- all_results_from_epi_stroma |> dplyr::filter(Zone == "Neck")
epi_stroma_neck 

# creating df of only differntially expressed genes in neck
epi_stroma_neck_sig <- epi_stroma_neck |> dplyr::filter(delabel != "NA")
epi_stroma_neck_sig

#vulcano plot isthmus only 

ggplot(data=epi_stroma_neck, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "log2 fold change", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


#create excel file 
write_excel_csv(epi_stroma_neck_sig, file = "epi_stroma_neck_sig.csv", ",")

# pathway analysis with ClusterProfiler for NECK ONLY

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")

colnames(epi_stroma_neck_sig) <- c("Gene", "Zone", "Contrast", "log2FoldChange", "pvalue", "FDR", "diffexpressed", "delabel")
epi_stroma_neck_sig

# log2 fold change
epi_stroma_neck_genelist <- epi_stroma_neck_sig$log2FoldChange

# name the vector
names(epi_stroma_neck_genelist) <- epi_stroma_neck_sig$Gene

# sort the list in decreasing order (required for clusterProfiler)
epi_stroma_neck_genelist = sort(epi_stroma_neck_genelist, decreasing = TRUE)
epi_stroma_neck_genelist

gse_epistr_neck <- gseGO(epi_stroma_neck_genelist, 
                        ont ="ALL", 
                        keyType = "SYMBOL", 
                        nPerm = 10000, 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = organism, 
                        pAdjustMethod = "none")

results_GSA_epistr_neck <- gse_epistr_neck@result
write_excel_csv(results_GSA_epistr_neck, file = "results_GSA_epistr_neck.csv", ",")
gse_epistr_neck@result

# Dotplot
require(DOSE)
dotplot(gse_epistr_neck, showCategory=30, split=".sign", label_format = 5) + facet_grid(.~.sign)

# Encrichment Map
emapplot(gse_epistr_neck, showCategory = 10) 
# Error: Error in has_pairsim(x) : 
# Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.

#install.packages("ggnewscale")

ema_neck <- pairwise_termsim(gse_epistr_neck, method = "JC", semData = NULL, showCategory = 200)
emapplot(ema_neck, showCategory = 10)

# Ridgeplot
ridgeplot(gse_epistr_neck) + labs(x = "enrichment distribution")

cnetplot(gse_epistr_neck, categorySize="pvalue", foldChange=epi_stroma_fov_genelist, showCategory = 3)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse_epistr_neck, by = "all", title = gse_epistr_neck$Description[3], geneSetID = 1)
 

# ============filtering and pathway analysis for BASE ====================

#filter for base

epi_stroma_base <- all_results_from_epi_stroma |> dplyr::filter(Zone == "Base")
epi_stroma_base


# creating df of only differntially expressed genes in base
epi_stroma_base_sig <- epi_stroma_base |> dplyr::filter(delabel != "NA")

#vulcano plot base only 

ggplot(data=epi_stroma_base, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "log2 fold change", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


#create excel file 

write_excel_csv(epi_stroma_base_sig, file = "epi_stroma_base_sig.csv", ",")

# pathway analysis with ClusterProfiler for BASE ONLY

colnames(epi_stroma_base_sig) <- c("Gene", "Zone", "Contrast", "log2FoldChange", "pvalue", "FDR", "diffexpressed", "delabel")
epi_stroma_base_sig

# log2 fold change
epi_stroma_base_genelist <- epi_stroma_base_sig$log2FoldChange

# name the vector
names(epi_stroma_base_genelist) <- epi_stroma_base_sig$Gene

# sort the list in decreasing order (required for clusterProfiler)
epi_stroma_base_genelist = sort(epi_stroma_base_genelist, decreasing = TRUE)
epi_stroma_base_genelist

gse_epistr_base <- gseGO(epi_stroma_base_genelist, 
                        ont ="ALL", 
                        keyType = "SYMBOL", 
                        nPerm = 10000, 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = organism, 
                        pAdjustMethod = "none")

results_GSA_epistr_base <- gse_epistr_base@result
write_excel_csv(results_GSA_epistr_base, file = "results_GSA_epistr_base.csv", ",")
gse_epistr_base@result

# Dotplot
require(DOSE)
dotplot(gse_epistr_base, showCategory=30, split=".sign", label_format = 5) + facet_grid(.~.sign)

# Encrichment Map
emapplot(gse_epistr_base, showCategory = 10) 
# Error: Error in has_pairsim(x) : 
# Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.

ema_base <- pairwise_termsim(gse_epistr_base, method = "JC", semData = NULL, showCategory = 200)
emapplot(ema_base, showCategory = 10)

# Ridgeplot

ridgeplot(gse_epistr_base) + labs(x = "enrichment distribution")

cnetplot(gse_epistr_base, categorySize="pvalue", foldChange=epi_stroma_fov_genelist, showCategory = 3)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse_epistr_base, by = "all", title = gse_epistr_base$Description[3], geneSetID = 1)

# =================================================================
# inter epithelium 

pData(target_nano_healthy_nQ3)$testZone <-
  factor(pData(target_nano_healthy_nQ3)$Zone, c("Foveola", "Isthmus", "Neck", "Base", "Muscularis"))

pData(target_nano_healthy_nQ3)[["patient"]] <-
  factor(pData(target_nano_healthy_nQ3)[["patient"]])

pData(target_nano_healthy_nQ3)[["class"]] <-
  factor(pData(target_nano_healthy_nQ3)[["class"]])

pData(target_nano_healthy_nQ3)[["segment"]] <-
  factor(pData(target_nano_healthy_nQ3)[["segment"]])

pData(target_nano_healthy_nQ3)[["Zone"]] <-
  factor(pData(target_nano_healthy_nQ3)[["Zone"]])

assayDataElement(object = target_nano_healthy_nQ3, elt = "log_q") <-
  assayDataApply(target_nano_healthy_nQ3, 2, FUN = log, base = 2, elt = "q_norm")

# =====
# Epithelium

ind <- pData(target_nano_healthy_nQ3)$class == "epithelium"
mixedOutmc <-
  mixedModelDE(target_nano_healthy_nQ3[, ind],
               elt = "log_q",
               modelFormula = ~ testZone + (1 + testZone | patient),
               groupVar = "testZone",
               nCores = parallel::detectCores(),
               multiCore = F)

# formatting the result to be more readable
r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(r_test)
r_test <- as.data.frame(r_test)
r_test$Contrast <- tests
r_test$Gene <- 
  unlist(lapply(colnames(mixedOutmc),
                rep, nrow(mixedOutmc["lsmeans", ][[1]])))
r_test$Class <- "epithelium" # note that I renamed this variable
r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
r_test <- r_test[, c("Gene", "Class", "Contrast", "Estimate", 
                     "Pr(>|t|)", "FDR")]

epithelium_comparisons <- r_test

#write table
write_excel_csv(epithelium_comparisons, file = "epithelium_comparisons.csv", ",")

epithelium_significant_comparisons <- r_test |> dplyr::filter(abs(Estimate) >= 0.6 & (Estimate) <= -0.6 & FDR <= 0.05)

# add a column of NAs
epithelium_comparisons$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
epithelium_comparisons$diffexpressed[epithelium_comparisons$Estimate > 0.6 & epithelium_comparisons$`Pr(>|t|)` < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN" 
epithelium_comparisons$diffexpressed[epithelium_comparisons$Estimate < -0.6 & epithelium_comparisons$`Pr(>|t|)` < 0.05] <- "DOWN"

# Now write down the name of genes beside the points.
# Create a new column "delabel" to my detaframe, that will contain the name of genes differentially expressed (NA in case they are not)
epithelium_comparisons$delabel <- NA
epithelium_comparisons$delabel[epithelium_comparisons$diffexpressed != "NO"] <- epithelium_comparisons$Gene[epithelium_comparisons$diffexpressed != "NO"]

# separate data for comparisons 

epi_fov_ist <- epithelium_comparisons |> dplyr::filter(Contrast == "Foveola - Isthmus")
epi_fov_ist

# loop to separate the data.frame 
unique_test <- unique(epithelium_comparisons$Contrast)
unique_test

for (i in unique_test){
  assign(paste0("df_", i), subset(epithelium_comparisons, Contrast == i))
  write_excel_csv(subset(epithelium_comparisons, Contrast == i), file = paste(i, ".xlsx"))
}

#vulcano plots 

#Foveola-base
ggplot(data=`df_Foveola - Base`, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "Base <-    -> Foveola", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#Foveola-isthmus
ggplot(data=`df_Foveola - Isthmus`, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "Isthmus <-    -> Foveola", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#Foveola-Neck
ggplot(data=`df_Foveola - Neck`, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "Neck <-    -> Foveola", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#Isthmus - Neck
ggplot(data=`df_Isthmus - Neck`, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "Neck <-    -> Isthmus", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#Isthmus - Base
ggplot(data=`df_Isthmus - Base`, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  labs(x = "Base <-    -> Isthmus", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#Neck - Base
ggplot(data=`df_Neck - Base`, aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  theme(legend.position="none")+
  geom_text_repel() +
  labs(x = "Base <-    -> Neck", y = "Significance, -log10(P-value)") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# is there a way to make a loop to create those graphs for all the conditions?
# loop to create graphs for all conditions

conditions_plots = list()
for(i in unique_test) {
  conditions_plots = ggplot((subset(epithelium_comparisons, Contrast == i)), aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=diffexpressed, label=delabel), show.legend = FALSE) +
    geom_point() + 
    theme_minimal() +
    theme(legend.position="none")+
    geom_text_repel() +
    ggtitle(str_c("Comparison (A-B) ", i)) +
    labs(x = "B <-    -> A", y = "Significance, -log10(P-value)") +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  print(conditions_plots)
}

# to fix details about the naming but loop works. 

# creating df of only differentially expressed genes with a for loop

for (i in unique_test){
  assign(paste0("df_sign_", i), subset(epithelium_comparisons, Contrast == i) |> dplyr::filter(delabel != "NA")) 
  write_excel_csv(subset(epithelium_comparisons, Contrast == i), file = paste(i, ".xlsx"))
}

# now I want to try to do other clusterprofiler analysis with the loop. 
# create a list of data.frames
intra_epithelium = list()
for (i in unique_test){
  intra_epithelium[[i]] <- subset(epithelium_comparisons, Contrast == i)
}

# create the filtered version 
intra_epithelium_filtr = list()
for (i in unique_test){
  intra_epithelium_filtr[[i]] <- subset(epithelium_comparisons, Contrast == i) |> dplyr::filter(delabel != "NA")
}

# for each df_sign_ 
# make a vector from the column Pr(>|t|)
intra_epi_genelists = list()
for (i in unique_test){
  assign(paste('vec_', i), intra_epithelium_filtr[[i]])
}

#from here.

# log2 fold change
epi_stroma_fov_genelist <- epi_stroma_foveola_sig$log2FoldChange

for (i in unique_test){
  assign(paste0("genlist_", i), subset(epithelium_comparisons, Contrast == i) |> dplyr::filter(delabel != "NA"))
}

# name the vector
names(epi_stroma_ist_genelist) <- epi_stroma_isthmus_sig$Gene

# sort the list in decreasing order (required for clusterProfiler)
epi_stroma_ist_genelist = sort(epi_stroma_ist_genelist, decreasing = TRUE)
epi_stroma_ist_genelist

gse_epistr_ist <- gseGO(epi_stroma_ist_genelist, 
                        ont ="ALL", 
                        keyType = "SYMBOL", 
                        nPerm = 10000, 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = organism, 
                        pAdjustMethod = "none")

results_GSA_epistr_ist <- gse_epistr_ist@result
write_excel_csv(results_GSA_epistr_ist, file = "results_GSA_epistr_ist.csv", ",")
gse_epistr_ist@result

# Dotplot
require(DOSE)
dotplot(gse_epistr_ist, showCategory=30, split=".sign", label_format = 5) + facet_grid(.~.sign)

# Encrichment Map
emapplot(gse_epistr_ist, showCategory = 10) 
# Error: Error in has_pairsim(x) : 
# Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.

#install.packages("ggnewscale")
library(ggnewscale)
ema_ist <- pairwise_termsim(gse_epistr_ist, method = "JC", semData = NULL, showCategory = 200)
emapplot(ema_ist, showCategory = 10)

# Ridgeplot
#install.packages("ggridges")
library(ggridges)
ridgeplot(gse_epistr_ist) + labs(x = "enrichment distribution")

cnetplot(gse_epistr_ist, categorySize="pvalue", foldChange=epi_stroma_fov_genelist, showCategory = 3)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse_epistr_ist, by = "all", title = gse_epistr_ist$Description[3], geneSetID = 1)



# =====
# Stroma

ind <- pData(target_nano_healthy_nQ3)$class == "stroma"
mixedOutmc <-
  mixedModelDE(target_nano_healthy_nQ3[, ind],
               elt = "log_q",
               modelFormula = ~ testZone + (1 + testZone | patient),
               groupVar = "testZone",
               nCores = parallel::detectCores(),
               multiCore = F)

# formatting the result to be more readable
r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(r_test)
r_test <- as.data.frame(r_test)
r_test$Contrast <- tests
r_test$Gene <- 
  unlist(lapply(colnames(mixedOutmc),
                rep, nrow(mixedOutmc["lsmeans", ][[1]])))
r_test$Class <- "stroma" # note that I renamed this variable
r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
r_test <- r_test[, c("Gene", "Class", "Contrast", "Estimate", 
                     "Pr(>|t|)", "FDR")]

stroma_comparisons <- r_test

#write table
write_excel_csv(stroma_comparisons, file = "stroma_comparisons.csv", ",")

stroma_significant_comparisons <- r_test |> dplyr::filter(abs(Estimate) >= 0.6 & (Estimate) <= -0.6 & FDR <= 0.05)




# Using biobroom
# Installing them first:
# if (!requireNamespace("BiocManager", quietly=TRUE))
# install.packages("BiocManager")
# BiocManager::install("biobroom")

