####################################################################
# Author    : M. Dubbelaar
# Date      : 15-sept-2015
# File Name : EdgeRVisualisation.r
# Purpose   : This script creates a visualisation of the data with 
#             the use of the EdgeR package.
# Used Files: EdgeRFunctions.R
#             plotColors.R
#             loadingEnsemblData.R
#             EdgeRLinearTime.R
#             EdgeRDispersions.R
#             EdgeRScatterplots.R
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
#install.packages("rgl")

library(limma)
library(edgeR)
library(biomaRt)
library(gplots)
library(rgl)
resultPathway <- "/home/mdubbelaar/Desktop/APP23_results/EdgeR/"
####################################################################
#               Reading of the data and the targets                #
####################################################################
source("Code/EdgeRFunctions.R")
M1 <- getData("/media/mdubbelaar/6CEC0BDEEC0BA186/1507_Holtman_RNAseq/run01/results/expression/expressionTable/expression_table02.genelevel.GRCm38.v76.htseq.txt.table")
targets <- getTarget("/home/mdubbelaar/APP23/Targets.csv")
source("../plotColors.R")

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
#                          APP Data Check                          #
####################################################################
# The CPMatrix contains the cpm values obtained from the dge vector.
CPMmatrix <- as.matrix(cpm(dge, log = T))
rownames(CPMmatrix) <- BioM[,3]

# The App gene is filtered from the CPMatrix. This is the gene that
# should show a difference when looking at the WT and HET mice.
appGene <- which(rownames(CPMmatrix)=="App")
# The values of this genes are obtained from the CPMatrix and saved
# as a plot that shows the expression of the App gene.
appExpression <- CPMmatrix[appGene,]

pdf(paste(resultPathway, "Plots/AppCheck.pdf", sep=""))
plot(appExpression, col=col_cell_age, pch=20, cex=3, ylab = "Expression of App")
dev.off()
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
# The M3 dataset is made so it can be used in the creation of the heatmaps
# Without interfering with the M2 dataset. The M3 dataset contains
# genenames instead of the Ensemble id's.
M3 <- M2
rownames(M3) <- BioM[,3]
####################################################################
#                            MDS Plot                              #
####################################################################
# The MDS plot is made and saved within a file with the use of pdf()
# dev.off() is necessary to close the connection to the document.
pdf(paste(resultPathway, "Plots/MDS_plot_age_effect.pdf", sep=""))
plotMDS.DGEList(dge, col=col_cell_age, main="MDS plot (age effect)")
dev.off() 
####################################################################
#                    Differential Expression                       #
####################################################################
source("Code/EdgeRLinearTime.R")
source("Code/EdgeRDispersions.R")
source("Code/EdgeRScatterplots.R")
setwd("/home/mdubbelaar/Desktop/APP23_results/EdgeR/")
source("/home/mdubbelaar/APP23/CreatedData/saveDEData.R")
setwd("/home/mdubbelaar/APP23/CreatedData/EdgeR/")
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
pdf(paste(resultPathway, "Plots/Spearman_Cor_Plot.pdf", sep=""))
heatmap.2(cor, symm = TRUE, col=greenred(350),
          trace="none", cexRow = 1 , cexCol = 1, ColSideColors= col_cell_age, RowSideColors=col_cell_age)
dev.off()
####################################################################
#                           PCA Plot                               #
####################################################################
pcaPlot(resultPathway)
####################################################################
#                            Heatmaps                              #
####################################################################
# Checks the two toptable results with the most unique genes.
# These heatmaps are saved within the heatmap.pdf.
pdf(paste(resultPathway, "Plots/heatmaps_containing_most_genes.pdf", sep=""))
heatmap.2(M2[match(rownames(toptable11.results), rownames(M2)),c(4:6, 10:12, 16:18, 22:24)], ColSideColors = col_cell_age[c(4:6, 10:12, 16:18, 22:24)], cexRow = 0.01, trace = "none", scale = "row", main="24M HET - 2M HET")
heatmap.2(M2[match(rownames(toptable12.results), rownames(M2)),c(1:3, 7:9, 13:15, 19:21)], ColSideColors = col_cell_age[c(1:3, 7:9, 13:15, 19:21)], cexRow = 0.01, trace = "none", scale = "row", main="24M WT - 2M WT")
dev.off()
