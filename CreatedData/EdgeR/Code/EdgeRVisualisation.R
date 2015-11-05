####################################################################
# Author    : M. Dubbelaar
# Date      : 15-sept-2015
# File Name : EdgeRVisualisation.r
# Purpose   : This script creates a visualisation of the data with 
#             the use of the EdgeR package.
# Used Files: loadingAppFile.R
#             loadingEnsemblData.R
#             plotColors.R
#             EdgeRLinearTime.R
#             EdgeRDispersions.R
#             EdgeRPCA.R
#             saveDEData.R
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/APP23/CreatedData/EdgeR/")
#biocLite("limma")
#biocLite("edgeR")
#biocLite("biomaRt")
#biocLite("gplots")
library(limma)
library(edgeR)
library(biomaRt)
library(gplots)
####################################################################
#               Reading of the data and the targets                #
####################################################################
source("../loadingAppFile.R")
# The data found in loadingAppFile needs some adjustment for further
# progress. The data needs to be filtered of the genenames within 
# the column and needs to be set as the column name.
M1 <- APP23_data[,2:25]
row.names(M1) <- APP23_data[,1]

# A dge list is creates, this list is used to estimate the GLM 
# common, trended and tagwise dispersion. This object holds the 
# dataset that needs to be analysed by EgdeR and the calculations 
# which where performed on the dataset.
dge <- DGEList(counts=M1, group=factor(targets$Conditie) )
####################################################################
#                      Loading Ensembl data                        #
####################################################################
source("../loadingEnsemblData.R")
####################################################################
#                        Necessary functions                       #
####################################################################
# This function is made to calculate the unique genes for the differential expressions.
# The genes with a FDR below a number of 0.05 are returned.
filterGenes <- function(tabel) {
  uniqueGenes <- which(tabel[[1]]$FDR < 0.05)
  uniqueGenes <- tabel[uniqueGenes,]
  return(uniqueGenes)
}
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
####################################################################
#                            MDS Plot                              #
####################################################################
source("../plotColors.R")
# The MDS plot is made and saved within a file with the use of pdf()
# dev.off() is necessary to close the connection to the document.
pdf("Plots/MDS_plot_age_effect.pdf") 
plotMDS.DGEList(dge, col=col_cell_age, main="MDS plot (age effect)")
dev.off() 
####################################################################
#                    Differential Expression                       #
####################################################################
source("Code/EdgeRLinearTime.R")
source("Code/EdgeRDispersions.R")
source("../saveDEData.R")
####################################################################
#                          Quality plots                           #
####################################################################
# A modelling of the mean-varance relationship in the DGE data 
# is made with the use of plotMeanVar. The gene means and variances
# are tagged and are shown within the plot.
# The show.raw.vars shows the found vars (the grey circles) and NBline
# shows the trend of the vars (the blue line).
plotMeanVar(dge, show.raw.vars = T, NBline = T)
# Plots the biological coefficient of variantion agains the log2 
# counts per million. The common, trended ans tagwise BCV estimates
# are shown in the plot
plotBCV(dge, cex=0.4)
# Plots the sample relations. The distances of the RNA-seq is 
# calulated, the distances represent the coefficient of variation
# of expression among the sampels.
plotMDS.DGEList(dge, col=col_cell_age)
####################################################################
#                    Spearman correlatie plot                      #
####################################################################
# A heatmap with the correlation of the different genes are shown.
# This correlation is calculated from the M2 dataset.
cor <- cor(M2, method="spearman")
rownames(cor) <- targets$Samples
colnames(cor) <- targets$Conditie
pdf("Plots/Spearman_Cor_Plot.pdf")
heatmap.2(cor, symm = TRUE, col=greenred(350),
          trace="none", cexRow = 1 , cexCol = 1, ColSideColors= col_cell_age, RowSideColors=col_cell_age)
dev.off()
####################################################################
#                           PCA Plot                               #
####################################################################
source("Code/EdgeRPCA.R")
####################################################################
#                            Heatmaps                              #
####################################################################
# Checks the two toptable results with the most unique genes.
# These heatmaps are saved within the heatmap.pdf.
pdf("Plots/heatmaps_containing_most_genes.pdf") 
heatmap.2(M2[match(rownames(toptable11.results), rownames(M2)),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="24M HET - 6-8W HET")
heatmap.2(M2[match(rownames(toptable12.results), rownames(M2)),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="24M WT - 6-8W WT")
dev.off()
