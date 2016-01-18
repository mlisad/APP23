####################################################################
# Author    : M. Dubbelaar
# Date      : 08-jan-2016
# File Name : EdgeRHeatmapsMerge.R
# Purpose   : Merging heatmaps from the GLM to create one heatmap
#             that can be used for the manuscript.
####################################################################
#                         Loading data (old)                       #
####################################################################
# The genes and the values are saved with the use of the function plotMostExpr.
dataMAO <- read.csv(paste(resultPathway, "Made_Documents/Top20Genes/ResultsForOneHeatmapTop20oldMA.txt", sep = ""))
dataMGO <- read.csv(paste(resultPathway, "Made_Documents/Top20Genes/ResultsForOneHeatmapTop20oldMG.txt", sep = ""))
dataIEO <- read.csv(paste(resultPathway, "Made_Documents/Top20Genes/ResultsForOneHeatmapTop20oldIE.txt", sep = ""))
####################################################################
dataMat <- getDataMat(dataMAO, dataMGO, dataIEO)
# The following colors are made to distinct the different effect.
# This gives a good overview were a specific gene is located from.
col_genes <- rep("black", length(dataMat[,1]))
col_genes[dataMat[,2] != 0] = col="green"
col_genes[dataMat[,2] != 0] = col="green"
col_genes[dataMat[,3] != 0] = col="purple"
col_genes[dataMat[,4] != 0] = col="grey"
col_genes[dataMat[,2] != 0 & dataMat[,4] != 0] = col = "palegreen"
col_genes[dataMat[,3] != 0 & dataMat[,4] != 0] = col = "orchid"
col_genes[dataMat[,2] != 0 & dataMat[,3] != 0] = col = "grey47"
col_genes[dataMat[,2] != 0 & dataMat[,3] != 0 & dataMat[,4] != 0] = col = "black"
# The color of the heatmap itself gets the colors; purple, lightgrey and green.
hmcol <- colorRampPalette(c("purple", "lightgrey", "green")) (n=99)
rownamesData <- dataMat[,1]
####################################################################
#                         LogFC heatmap                            #
####################################################################
logFCVal <- getlogFCVal(rownamesData, 7:24, 18, 19)
# The colors for the genes are cut off into another matrix.
logFCResults <- logFCVal[,1:18]
colnames(logFCResults) <- targets$Conditie[7:24]
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
pdf(paste(resultPathway, "Plots/HeatmapMergedOld.pdf", sep=""))
heatmap.2(logFCResults, RowSideColors = logFCVal[,19], ColSideColors= col_side, col=hmcol, trace="none", dendrogram = "none",Colv = F, Rowv = F, cexRow = 0.7, labCol = F, density.info="none")
dev.off()


####################################################################
#                       Loading data (young)                       #
####################################################################
dataMAY <- read.csv(paste(resultPathway, "Made_Documents/Top20Genes/ResultsForOneHeatmapTop20youngMA.txt", sep = ""))
dataMGY <- read.csv(paste(resultPathway, "Made_Documents/Top20Genes/ResultsForOneHeatmapTop20youngMG.txt", sep = ""))
dataIEY <- read.csv(paste(resultPathway, "Made_Documents/Top20Genes/ResultsForOneHeatmapTop20youngIE.txt", sep = ""))
####################################################################
dataMat <- getDataMat(dataMAY, dataMGY, dataIEY)
# The following colors are made to distinct the different effect.
# This gives a good overview were a specific gene is located from.
col_genes <- rep("black", length(dataMat[,1]))
col_genes[dataMat[,2] != 0] = col="green"
col_genes[dataMat[,2] != 0] = col="green"
col_genes[dataMat[,3] != 0] = col="purple"
col_genes[dataMat[,4] != 0] = col="grey"
col_genes[dataMat[,2] != 0 & dataMat[,4] != 0] = col = "palegreen"
col_genes[dataMat[,3] != 0 & dataMat[,4] != 0] = col = "orchid"
col_genes[dataMat[,2] != 0 & dataMat[,3] != 0] = col = "grey47"
col_genes[dataMat[,2] != 0 & dataMat[,3] != 0 & dataMat[,4] != 0] = col = "black"
# The color of the heatmap itself gets the colors; purple, lightgrey and green.
hmcol <- colorRampPalette(c("purple", "lightgrey", "green")) (n=99)
rownamesData <- dataMat[,1]
####################################################################
#                         LogFC heatmap                            #
####################################################################
logFCVal <- getlogFCVal(rownamesData, 1:12, 12, 13)
# The colors for the genes are cut off into another matrix.
logFCResults <- logFCVal[,1:12]
colnames(logFCResults) <- targets$Conditie[1:12]
# The class needs to be resetted to numeric otherwise it won't create a heatmap.
class(logFCResults) <- "numeric"
####################################################################
# The colors of the column are made to distingious the genotype of the mice.
col_side <- rep("black", length(colnames(logFCResults)))
col_side[colnames(logFCResults) ==  "WT.02"] = col="lightpink"
col_side[colnames(logFCResults) ==  "WT.06"] = col="indianred1"
col_side[colnames(logFCResults) ==  "HET.02"] = col="cyan2"
col_side[colnames(logFCResults) ==  "HET.06"] = col="deepskyblue2"
####################################################################
# Creation of the heatmap as a pdf file.
pdf(paste(resultPathway, "Plots/HeatmapMergedYoung.pdf", sep=""))
heatmap.2(logFCResults, RowSideColors = logFCVal[,13], ColSideColors= col_side, col=hmcol, trace="none", dendrogram = "none",Colv = F, Rowv = F, cexRow = 0.7, labCol = F, density.info="none")
dev.off()