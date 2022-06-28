# install.packages("devtools")
# devtools::install_github("Nanostring-Biostats/NanoStringNCTools")
# devtools::install_github("Nanostring-Biostats/GeomxTools", ref = "dev")
# devtools::install_github("Nanostring-Biostats/GeoMxWorkflows", ref = "main")

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(readxl)

setwd("~/OneDrive - Norwich BioScience Institutes/Documents/nanostring/Nanostring_CZ_MO/")
datadir <- file.path("~/OneDrive - Norwich BioScience Institutes/Documents/nanostring/Nanostring_CZ_MO")
DCCdir <- file.path("~/OneDrive - Norwich BioScience Institutes/Documents/nanostring/Nanostring_CZ_MO/GeoMxNGS_Pipeline_Requeue_Normal_samples_CZ_DCC")

DCCFiles <- dir(file.path(datadir, "GeoMxNGS_Pipeline_Requeue_Normal_samples_CZ_DCC"),
  pattern = ".dcc$",
  full.names = TRUE, recursive = TRUE
)

PKCFiles <- file.path(datadir, "Hs_R_NGS_WTA_v1.0.pkc")

SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

nano_healthy <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = SampleAnnotationFile,
  phenoDataSheet = "Annotation template",
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = c("aoi", "roi"),
  experimentDataColNames = c("panel")
)

# Asking which object is
class(nano_healthy)

isS4(nano_healthy)
is(nano_healthy, "ExpressionSet")

# access the count matrix
assayData(nano_healthy)[["exprs"]][1:3, 1:3] # genes and count per each DCC file
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

# Adding the pdata for the zones, using the same method

zone <- rep(c("Foveola", "Foveola", "Isthmus", "Isthmus", "Neck", "Neck", "Base", "Base", "Muscularis"), times = 9)

pData(nano_healthy)[, c("Zone")] <- zone
pData(nano_healthy)[, c("Zone")]

# Adding an extra column for the simplified names of the patients.

patient <- rep(c("27C", "29C", "31C"), each = 27)
patient
pData(nano_healthy)[, c("patient")] <- patient

pData(nano_healthy)

# Adding an extra column for epithelium vs stroma only.
# in this way we can do comparisons among the two classes only if we want to go easier
# we can also find clusters only upon those characteristics

class <- rep(c("epithelium", "stroma", "stroma", "epithelium", "stroma", "epithelium", "stroma", "epithelium", "stroma"), times = 9)
class
pData(nano_healthy)[, c("class")] <- class

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
count_mat <- count(pData(nano_healthy), `patient`, `class`, `segment`, `Zone`)

# gather the data and plot in order: class, slide name, region, segment
test_gr <- gather_set_data(count_mat, 1:4)
test_gr$x <- factor(test_gr$x,
  levels = c("patient", "class", "segment", "Zone")
)

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
nano_healthy_shifted <- shiftCountsOne(nano_healthy, useDALogic = TRUE) # you do not want to use zero values

# Segment QC
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below

# Why did i change those values

QC_params <-
  list(
    minSegmentReads = 1000, # Minimum number of reads (1000)
    percentTrimmed = 80, # Minimum % of reads trimmed (80%)
    percentStitched = 80, # Minimum % of reads stitched (80%)
    percentAligned = 80, # Minimum % of reads aligned (80%)
    percentSaturation = 50, # Minimum sequencing saturation (50%)
    minNegativeCount = 1, # Minimum negative control counts (10) #this because we have low negative counts
    maxNTCCount = 1000, # Maximum counts observed in NTC well (1000)
    minNuclei = 20, # Minimum # of nuclei estimated (100) #in vimentin we have low nuclei. Nanostring specialst suggested to keep anything >20
    minArea = 1000
  ) # Minimum segment area (5000)


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
# QC_histogram(nano_healthy_shifted, "nuclei", col_by, 20)

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
# there are different fragments that fall behind the dotted line. However, I am still ot sure about the meaning of the dotted line.

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

# exclude oulier probes

# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <-
  subset(
    nano_healthy_shift_QC,
    fData(nano_healthy_shift_QC)[["QCFlags"]][, c("LowProbeRatio")] == FALSE &
      fData(nano_healthy_shift_QC)[["QCFlags"]][, c("GlobalGrubbsOutlier")] == FALSE
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
for (module in modules) {
  vars <- paste0(
    c("NegGeoMean_", "NegGeoSD_"),
    module
  )
  if (all(vars[1:2] %in% colnames(pData(target_nano_healthy_SQC)))) {
    LOQ[, module] <-
      pmax(
        minLOQ,
        pData(target_nano_healthy_SQC)[, vars[1]] *
          pData(target_nano_healthy_SQC)[, vars[2]]^cutoff
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
    limits = c(0, 1),
    labels = scales::percent
  ) +
  theme_bw() +
  scale_y_continuous(
    labels = scales::percent, limits = c(0, 1),
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

library(reshape2) # for melt
library(cowplot) # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives

exprs(target_nano_healthy_norm)

ann_of_interest <- "class"
stat_norm <-
  data.frame(
    row.names = colnames(exprs(target_nano_healthy_norm)),
    Segment = colnames(exprs(target_nano_healthy_norm)),
    Annotation = pData(nano_healthy)[, ann_of_interest],
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
  rel_widths = c(0.43, 0.57)
)
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

#### NORMALISATION####
# I can create two dataset, one normalized Q3 and another one for background

#### Q3 norm (75th percentile) for WTA  with or without custom spike-ins
target_nano_healthy_nQ3 <-NanoStringNCTools::normalize(target_nano_healthy_norm,
  #data_type = "RNA",
  norm_method = "quant",
  desiredQuantile = .75,
  toElt = "q_norm"
)

# Background normalization for WTA/CTA without custom spike-in
target_nano_healthy_backg <- normalize(target_nano_healthy_norm,
#  data_type = "RNA",
  norm_method = "neg",
  fromElt = "exprs",
  toElt = "neg_norm"
)

# visualize the first 15 segments with each normalization method
boxplot(exprs(target_nano_healthy_nQ3)[, 1:15],
  col = "#9EDAE5", main = "Raw Counts",
  log = "y", names = 1:15, xlab = "Segment",
  ylab = "Counts, Raw"
)


boxplot(assayDataElement(target_nano_healthy_nQ3[, 1:15], elt = "q_norm"),
  col = "#2CA02C", main = "Q3 Norm Counts",
  log = "y", names = 1:15, xlab = "Segment",
  ylab = "Counts, Q3 Normalized"
)

boxplot(assayDataElement(target_nano_healthy_backg[, 1:15], elt = "neg_norm"),
  col = "#FF7F0E", main = "Neg Norm Counts",
  log = "y", names = 1:15, xlab = "Segment",
  ylab = "Counts, Neg. Normalized"
)

# the plots are slightly different but I am not sure which one is the best to use.


### unsupervised analysis
#install.packages("umap")
library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
# what does happen if I add or remove seeds?
# Marton: You will get a slightly different output every time. Selecting a seed will
# allow you to be consistent within a run at least.
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

umap_out <-
  umap(t(log2(assayDataElement(target_nano_healthy_nQ3, elt = "q_norm"))),
    config = custom_umap
  )
pData(target_nano_healthy_nQ3)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1, 2)]

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

# change seeds to 50 to see what changes

custom_umap1 <- umap::umap.defaults
custom_umap1$random_state <- 50

umap_out1 <-
  umap(t(log2(assayDataElement(target_nano_healthy_nQ3, elt = "q_norm"))),
    config = custom_umap1
  )
pData(target_nano_healthy_nQ3)[, c("UMAP1", "UMAP2")] <- umap_out1$layout[, c(1, 2)]
ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = UMAP1, y = UMAP2, color = Zone, shape = segment)
) +
  geom_point(size = 3) +
  theme_bw() # graph changes completely.
# need to understand which are the best conditions.

# try with 3 components and 42 seeds

custom_umap2 <- umap::umap.defaults
custom_umap2$random_state <- 42
custom_umap2$n_components <- 3

umap_out2 <-
  umap(t(log2(assayDataElement(target_nano_healthy_nQ3, elt = "q_norm"))),
    config = custom_umap2
  )
pData(target_nano_healthy_nQ3)[, c("UMAP1", "UMAP2", "UMAP3")] <- umap_out2$layout[, c(1, 2, 3)]
ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = UMAP1, y = UMAP2, z = UMAP3, color = Zone, shape = segment)
) +
  geom_point(size = 3) +
  theme_bw()

# Do not know how to add the z axes to the plot #however with the same number of seeds I got a different type of graph, why?

# now with the background normalisation

umap_out_BG <-
  umap(t(log2(assayDataElement(target_nano_healthy_backg, elt = "neg_norm"))),
    config = custom_umap
  )
pData(target_nano_healthy_backg)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1, 2)]

Bkg_umap <- ggplot(
  pData(target_nano_healthy_backg),
  aes(x = UMAP1, y = UMAP2, color = Zone, shape = segment)
) +
  geom_point(size = 3) +
  theme_bw()


# plot the umap from the two normalisations one next to the other to see the difference

library(ggpubr)

ggarrange(Bkg_umap + rremove("legend"), Q3_umap, labels = c("background norm", "Q3 norm"), ncol = 2, nrow = 1)
# to improve the position of the lables, and the distribution of the graphs.

# it seems that the two normalisation methods do not give particular differences in the segregation of the samples.

# Marton: correct!

# run tSNE
# Marton: not sure tSNE is needed if you already did UMAP
# seed 42, on Q3 normalised data and background.
# plotting segments and zones

set.seed(42) # set the seed for tSNE as well

# Q3 normalisation
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_nano_healthy_nQ3, elt = "q_norm"))),
    perplexity = ncol(target_nano_healthy_nQ3) * .15
  )
pData(target_nano_healthy_nQ3)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1, 2)]

Q3_tsne <- ggplot(
  pData(target_nano_healthy_nQ3),
  aes(x = tSNE1, y = tSNE2, color = segment, shape = zone)
) +
  geom_point(size = 3) +
  theme_bw()


# same for the normalised data upon background noise

tsne_out1 <-
  Rtsne(t(log2(assayDataElement(target_nano_healthy_backg, elt = "neg_norm"))),
    perplexity = ncol(target_nano_healthy_backg) * .15
  )
pData(target_nano_healthy_backg)[, c("tSNE1", "tSNE2")] <- tsne_out1$Y[, c(1, 2)]

Bkg_tsne <- ggplot(
  pData(target_nano_healthy_backg),
  aes(x = tSNE1, y = tSNE2, color = segment, shape = zone)
) +
  geom_point(size = 3) +
  theme_bw()

# putting all the graph together.

ggarrange(Bkg_umap + rremove("legend"), Q3_umap, Bkg_tsne + rremove("legend"), Q3_tsne,
  labels = c("background norm", "Q3 norm"), ncol = 2, nrow = 2
)

# to be improved: removing the legend and putting it as an extra so that all the graphs will appear equal.
# the two graphs have a differnet color annotation, so I cannot use the same legend for all the graphs.

#### CLUSTERING HIGH CV GENES####

# from here we might not need it since we might use the counts for our network analysis

library(pheatmap) # for heatmap
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
  clustering_method = "complete", # I tryied to change the method but it seems that the clusters do not make much sense.
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  breaks = seq(-3, 3, 0.05),
  color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
  annotation_col =
    pData(target_nano_healthy_nQ3)[, c("Zone", "segment", "patient", "class")]
) # how can I change the colours of the annotation?
#Marton: I changed it to complete linkage.
# Complete linkage clustering means that the distance from one cluster to another is calculated based on the furthest members of the cluster. The used clustering is sensitive for the furthest elements. Complete linkage does not join together with the furthest clusters, producing a clear picture.

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

# I am not sure I wrote the formula correctly. I need to understand the terms I need to insert in the formula.

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

# Keep significant only

all_results_from_epi_stroma_significant <- all_results_from_epi_stroma |> dplyr::filter(abs(Estimate) >= 1 & FDR <= 0.05)

# Results are all upregulated - the downregulated genes are not significant
# Plot below show significantly upregulated genes

ggplot(all_results_from_epi_stroma_significant, 
       aes(x = Gene, y = Estimate, fill=Zone)) +
  geom_bar(stat='identity')+
  coord_flip()+
  facet_wrap(~Zone)

# Violin plot for example gene (GZMB)

ggplot(pData(target_nano_healthy_nQ3),
       aes(x = class, fill = class,
           y = assayDataElement(target_nano_healthy_nQ3["GZMB", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "GZMB Expression") +
  scale_y_continuous(trans = "log2") +
  facet_wrap(~Zone) +
  theme_bw()

# between slides?
# Marton: To do the comparison between healthy & diseased we need to load the cancer data 
# Marton: If I understand correctly so far we have only loaded the healthy data

