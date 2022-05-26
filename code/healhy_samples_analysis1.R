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


