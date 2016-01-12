####################################################################
# Author    : M. Dubbelaar
# Date      : 08-jan-2016
# File Name : EdgeRHeatmapsMerge.R
# Purpose   : Merging heatmaps from the GLM to create one heatmap
#             that can be used for the manuscript.
####################################################################
#                           Loading data                           #
####################################################################
# The genes and the values are saved with the use of the function plotMostExpr.
dataMAO <- read.csv("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/ResultsForOneHeatmapOldMA.txt")
dataMGO <- read.csv("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/ResultsForOneHeatmapOldMG.txt")
dataIEO <- read.csv("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/ResultsForOneHeatmapOldIE.txt")

# The column names contain the condition of the samples
colnames(M3) <- targets$Conditie
# The unique genes are loaded from the three files are saved are rownamesData.
rownamesData <- unique(c(as.matrix(dataMGO[,1]), as.matrix(dataMAO[,1]), as.matrix(dataIEO[,1])))
# The values are loaded into a matrix, this matrix is made to detect from which group
# the gene originate from.
dataMat <- apply(cbind(rownamesData, dataMGO[match(rownamesData, dataMGO[,1]),2], dataMAO[match(rownamesData, dataMAO[,1]),2], dataIEO[match(rownamesData, dataIEO[,1]),2]), 2, function(x) unname(unlist(x)))
# The information containing NA is overwritten as 0.
dataMat[is.na(dataMat)] <- 0

####################################################################
# The following colors are made to distinct the different effect.
# This gives a good overview were a specific gene is located from.
col_genes <- rep("black", length(dataMat[,1]))
col_genes[dataMat[,2] != 0] = col="dodgerblue3"
col_genes[dataMat[,3] != 0] = col="green4"
col_genes[dataMat[,4] != 0] = col="orchid4"
col_genes[dataMat[,2] != 0 & dataMat[,4] != 0] = col = "springgreen"
col_genes[dataMat[,3] != 0 & dataMat[,4] != 0] = col = "lightcoral"
####################################################################
# The color of the heatmap itself gets the colors; purple, lightgrey and green.
hmcol <- colorRampPalette(c("purple", "lightgrey", "green")) (n=99)

# The real data is obtained with the use of the M3 dataset
# The values of the obtained data will be created with the use of the 
# raw information - the mean of the row. This creates the logFC.
logFCVal <- M3[match(rownamesData, rownames(M3)),7:24] - apply(M3[match(rownamesData, rownames(M3)),7:24], 1, mean)
# The colors for the genes are also added into the matrix
logFCVal <- apply(cbind(logFCVal, col_genes), 2, function(x) unname(unlist(x)))
# The rownames will be the genes.
rownames(logFCVal) <- rownames(M3[match(rownamesData, rownames(M3)),7:24] )
# An ordering is done on the group were the genes come from and the highest 
# expression of the HET24 mice.
logFCVal <- logFCVal[order(logFCVal[,19], logFCVal[,18], decreasing = T), order(nchar(colnames(logFCVal)))]
# The colors for the genes are cut off into another matrix.
logFCResults <- logFCVal[,1:18]
# The class needs to be resetted to numeric otherwise it won't create a heatmap.
class(logFCResults) <- "numeric"
####################################################################
# The colors of the column are made to distingious the genotype of the mice.
col_side <- rep("black", length(colnames(logFCResults)))
col_side[colnames(logFCResults) ==  "WT.06"] = col="indianred1"
col_side[colnames(logFCResults) ==  "WT.18"] = col="red"
col_side[colnames(logFCResults) ==  "WT.24"] = col="brown4"
col_side[colnames(logFCResults) ==  "HET.06"] = col="deepskyblue2"
col_side[colnames(logFCResults) ==  "HET.18"] = col="dodgerblue3"
col_side[colnames(logFCResults) ==  "HET.24"] = col="blue3"
####################################################################
# Creation of the heatmap as a pdf file.
pdf(paste(resultPathway, "Plots/HeatmapMerged.pdf", sep=""))
heatmap.2(logFCResults, RowSideColors = logFCVal[,19], ColSideColors= col_side, col=hmcol, trace="none", dendrogram = "none",Colv = F, Rowv = F, cexRow = 0.8, labCol = F, density.info="none")
dev.off()