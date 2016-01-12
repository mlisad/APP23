####################################################################
# Author    : M. Dubbelaar
# Date      : 15-sept-2015
# File Name : EdgeRVisualisation.r
# Purpose   : This script creates a visualisation of the data with 
#             the use of the EdgeR package.
# Used Files: EdgeRFunctions.R
#             plotColorsAPP23.R
#             plotColorsCKp25.R
#             plotColorsHuman.R
#             loadingEnsemblData.R
#             EdgeRLinearTimeAPP23.R
#             EdgeRLinearTimeCKp25.R
#             EdgeRDispersionsAPP23.R
#             EdgeRDispersionsCkp25.R
#             EdgeRDispersionsHuman.R
#             EdgeRPCA.R
#             saveDEData.R
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/APP23/CreatedData/Test/")
#biocLite("limma")
#biocLite("edgeR")
#biocLite("biomaRt")
#biocLite("gplots")
#install.packages("rgl")
#install.packages("gplots")
library(rgl)
library(limma)
library(edgeR)
library(biomaRt)
library(gplots)
library(data.table)
####################################################################
#                             Functions                            #
####################################################################
source("EdgeRFunctions.R")
####################################################################
#M1 <- getData("/media/mdubbelaar/BD_4T/APP23/HTSeq/Human/mergedCounts.txt")
#targets <- getTarget("/media/mdubbelaar/BD_4T/APP23/HTSeq/Human/TargetHuman.csv")
#col_cell_age <- rep("black", dim(targets)[1])
#col_cell_age[targets$Genotype ==  "CTRL"] = col="red"
#col_cell_age[targets$Conditie ==  "LOAD"] = col="blue3"
####################################################################
#M1 <- getData("/media/mdubbelaar/6CEC0BDEEC0BA186/1507_Holtman_RNAseq/run01/results/expression/expressionTable/expression_table02.genelevel.GRCm38.v76.htseq.txt.table")
#targets <- getTarget("/home/mdubbelaar/APP23/Targets.csv")
#source("/home/mdubbelaar/APP23/CreatedData/plotColors.R")
####################################################################
M1 <- getData("/media/mdubbelaar/BD_4T/APP23/HTSeq/Mouse/mergedCounts.txt")
targets <- getTarget("/media/mdubbelaar/BD_4T/APP23/HTSeq/Mouse/TargetMouse.csv")
col_cell_age <- rep("black", dim(targets)[1])
col_cell_age[targets$Conditie ==  "WT_02w"] = col="lightpink"
col_cell_age[targets$Conditie ==  "WT_06w"] = col="red"
col_cell_age[targets$Conditie ==  "CKp25_02w"] = col="deepskyblue2"
col_cell_age[targets$Conditie ==  "CKp25_06w"] = col="blue3"
####################################################################
#                        Filtering the data                        #
####################################################################
# A dge list is creates, this list is used to estimate the GLM 
# common, trended and tagwise dispersion. This object holds the 
# dataset that needs to be analysed by EgdeR and the calculations 
# which where performed on the dataset
dge <- DGEList(counts=M1, group=factor(targets$Conditie) )
####################################################################
#                      Loading Ensembl data                        #
####################################################################
source("/home/mdubbelaar/APP23/CreatedData/loadingEnsemblData.R")
####################################################################
#                       Estimation Dispersions                     #
####################################################################
# The data is filtered with a cut-off of 1 count per million (cpm).
# The data with 2 replicated samples will be saves into the vector.
isExpr <- rowSums(cpm(dge)>1) >= 2
# These expressed genes will be saved within the dge and BioM dataset.
dge <- dge[isExpr, ]
BioM <- BioM[isExpr,]

# The raw library sizes are calculated with the use of calNormFactors
# This is a step before calculating the estimates.
dge <- calcNormFactors(dge)
colnames(dge$counts) <- targets$Samples
M2 <- cpm(dge, log=TRUE)
M3 <- M2
rownames(M3) <- BioM[,3]
####################################################################
#                      Differential expression                     #
####################################################################
#differentialExpression("Human", "/home/mdubbelaar/Desktop/APP23_results/Test/Human/")
#differentialExpression("APP23", "/home/mdubbelaar/Desktop/APP23_results/Test/APP23/")
differentialExpression("CKp25", "/home/mdubbelaar/Desktop/APP23_results/Test/CKp25/")
