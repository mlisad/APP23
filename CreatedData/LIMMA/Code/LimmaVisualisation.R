####################################################################
# Author    : M. Dubbelaar
# Date      : 25-sept-2015
# File Name : LimmaVisualisation.R
# Purpose   : To visualise the APP data with the package limma.
# Used Files: loadingAppFile.R
#             loadingEnsemblData.R
#             LimmaLinearTime.R
#             LimmaDispersions.R
#             plotColors.R
#             saveDEData.R
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/APP23/CreatedData/LIMMA/")
#install.packages("gplots")
#install.packages("statmod")
#biocLite("biomaRt")
#biocLite("limma")
#biocLite("edgeR")
library(biomaRt)
library(gplots)
library(limma)
library(edgeR)
library(statmod)
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
dge <- DGEList(counts=M1)
M2 <- cpm(dge, log=TRUE)
####################################################################
#                      Loading Ensembl data                        #
####################################################################
source("../loadingEnsemblData.R")
####################################################################
#                          Filtering Data                          #
####################################################################
# The data is filtered with a cut-off of 1 count per million (cpm).
# The data with 2 replicated samples will be saves into the vector.
isExpr <- rowSums(cpm(dge)>1) >= 2
dge <- dge[isExpr, keep.lib.sizes=F]
# These expressed genes will be saved within the dge and BioM dataset.
BioM <- BioM[isExpr,]
# The raw library sizes are calculated with the use of calNormFactors
# This is a step before calculating the estimates.
dge <- calcNormFactors(dge)
####################################################################
#                        Necessary functions                       #
####################################################################
# This function is made to calculate the unique genes for the differential expressions.
# The genes with a adjusted p value < 0.05 are returned.
filterGenesWithLimma <- function(tabel) {
  uniqueGenes <- which(tabel$adj.P.Val < 0.05)
  uniqueGenes <- tabel[uniqueGenes,]
  return(uniqueGenes)
}
####################################################################
#                             MDS Plot                             #
####################################################################
# The MDS plot is made and saved within a file with the use of pdf()
# dev.off() is necessary to close the connection to the document.
pdf("/home/mdubbelaar/Desktop/APP23_results/LIMMA/Plots/MDS_plot.pdf")
plotMDS(dge, labels = targets$Conditie, col=as.numeric(targets$Conditie))
legend("topright", legend = unique(targets$Conditie), col=unique(targets$Conditie), pch=15, cex=.6)
dev.off()
####################################################################
#                              Colors                              #
####################################################################
source("../plotColors.R")
####################################################################
#                      Differential Expression                     #
####################################################################
source("Code/LimmaLinearTime.R")
source("Code/LimmaDispersions.R")
source("../saveDEData.R")
####################################################################
#                            Heatmaps                              #
####################################################################
# Checks the two toptable results with the most unique genes.
# These heatmaps are saved within the heatmap.pdf.
pdf("/home/mdubbelaar/Desktop/APP23_results/LIMMA/Plots/heatmaps_contains_most_unique.pdf")
heatmap.2(M2[match(rownames(toptable7.results), rownames(M2)),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="24M HET-18M HET")
heatmap.2(M2[match(rownames(toptable11.results), rownames(M2)),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="24M HET- 6-8W HET")
dev.off()
