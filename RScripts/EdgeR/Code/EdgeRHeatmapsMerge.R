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
#dataMAO <- read.csv("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/ResultsForOneHeatmapOldMA.txt")
#dataMGO <- read.csv("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/ResultsForOneHeatmapOldMG.txt")
#dataIEO <- read.csv("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/ResultsForOneHeatmapOldIE.txt")
dataMAO <- read.csv("/home/mdubbelaar/Desktop/ResultsForOneHeatmapTop20oldMA.txt")
dataMGO <- read.csv("/home/mdubbelaar/Desktop/ResultsForOneHeatmapTop20oldMG.txt")
dataIEO <- read.csv("/home/mdubbelaar/Desktop/ResultsForOneHeatmapTop20oldIE.txt")

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
col_genes[dataMat[,2] != 0] = col="green"
col_genes[dataMat[,3] != 0] = col="palegreen"
col_genes[dataMat[,4] != 0] = col="grey"
col_genes[dataMat[,2] != 0 & dataMat[,4] != 0] = col = "orchid"
col_genes[dataMat[,3] != 0 & dataMat[,4] != 0] = col = "purple"
col_genes[dataMat[,2] != 0 & dataMat[,3] != 0 & dataMat[,4] != 0] = col = "black"
####################################################################
# The color of the heatmap itself gets the colors; purple, lightgrey and green.
hmcol <- colorRampPalette(c("purple", "lightgrey", "green")) (n=99)
####################################################################
#                         LogFC heatmap                            #
####################################################################
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
logFCVal <- logFCVal[order(as.numeric(as.character(logFCVal[,18])), decreasing = T),]
logFCVal <- logFCVal[order(factor(logFCVal[,19], levels = c("green", "palegreen", "grey", "orchid", "purple", "black"))),]
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
pdf(paste(resultPathway, "Plots/HeatmapMergedOld.pdf", sep=""))
heatmap.2(logFCResults, RowSideColors = logFCVal[,19], ColSideColors= col_side, col=hmcol, trace="none", dendrogram = "none",Colv = F, Rowv = F, cexRow = 0.7, labCol = F, density.info="none")
dev.off()


####################################################################
#                     Z-score Heatmap (old)                        #
#                       Work in progress                           #
####################################################################
hmcol <- colorRampPalette(c("purple", "lightgrey", "green")) (n=99)
# The real data is obtained with the use of the M3 dataset
# The values of the obtained data will be created with the use of the 
# raw information - the mean of the row. This creates the logFC.
GeneInfo <- M3[match(rownamesData, rownames(M3)),7:24]
View(GeneInfo)
# The colors for the genes are also added into the matrix
GeneInfo <- apply(cbind(GeneInfo, col_genes), 2, function(x) unname(unlist(x)))
# The rownames will be the genes.
rownames(GeneInfo) <- rownames(M3[match(rownamesData, rownames(M3)),7:24] )
# An ordering is done on the group were the genes come from and the highest 
# expression of the HET24 mice.
GeneInfo <- GeneInfo[order(as.numeric(as.character(GeneInfo[,18])), decreasing = T),]
GeneInfo <- GeneInfo[order(factor(GeneInfo[,19], levels = c("green", "palegreen", "grey", "orchid", "purple", "black"))),]
# The colors for the genes are cut off into another matrix.
GeneResults <- GeneInfo[,1:18]
# The class needs to be resetted to numeric otherwise it won't create a heatmap.
class(GeneResults) <- "numeric"
####################################################################
# The colors of the column are made to distingious the genotype of the mice.
col_side <- rep("black", length(colnames(GeneResults)))
col_side[colnames(GeneResults) ==  "WT.06"] = col="indianred1"
col_side[colnames(GeneResults) ==  "WT.18"] = col="red"
col_side[colnames(GeneResults) ==  "WT.24"] = col="brown4"
col_side[colnames(GeneResults) ==  "HET.06"] = col="deepskyblue2"
col_side[colnames(GeneResults) ==  "HET.18"] = col="dodgerblue3"
col_side[colnames(GeneResults) ==  "HET.24"] = col="blue3"
####################################################################
# Creation of the heatmap as a pdf file.
pdf(paste(resultPathway, "Plots/HeatmapMergedOld_ZScores.pdf", sep=""))

heatmap.2(GeneResults, RowSideColors = GeneInfo[,19], ColSideColors= col_side, col=hmcol, trace="none", dendrogram = "none", Colv = F, Rowv = F, cexRow = 0.7, labCol = F, density.info="none", scale = "row")
dev.off()

####################################################################
#                       Loading data (young)                       #
####################################################################
dataMAY <- read.csv("/home/mdubbelaar/Desktop/ResultsForOneHeatmapTop20youngMA.txt")
dataMGY <- read.csv("/home/mdubbelaar/Desktop/ResultsForOneHeatmapTop20youngMG.txt")
dataIEY <- read.csv("/home/mdubbelaar/Desktop/ResultsForOneHeatmapTop20youngIE.txt")

# The column names contain the condition of the samples
colnames(M3) <- targets$Conditie
# The unique genes are loaded from the three files are saved are rownamesData.
rownamesData <- unique(c(as.matrix(dataMGY[,1]), as.matrix(dataMAY[,1]), as.matrix(dataIEY[,1])))
# The values are loaded into a matrix, this matrix is made to detect from which group
# the gene originate from.
dataMat <- apply(cbind(rownamesData, dataMGY[match(rownamesData, dataMGY[,1]),2], dataMAY[match(rownamesData, dataMAY[,1]),2], dataIEY[match(rownamesData, dataIEY[,1]),2]), 2, function(x) unname(unlist(x)))
# The information containing NA is overwritten as 0.
dataMat[is.na(dataMat)] <- 0
####################################################################
# The following colors are made to distinct the different effect.
# This gives a good overview were a specific gene is located from.
col_genes <- rep("black", length(dataMat[,1]))
col_genes[dataMat[,2] != 0] = col="green"
col_genes[dataMat[,3] != 0] = col="palegreen"
col_genes[dataMat[,4] != 0] = col="grey"
col_genes[dataMat[,2] != 0 & dataMat[,4] != 0] = col = "orchid"
col_genes[dataMat[,3] != 0 & dataMat[,4] != 0] = col = "purple"
col_genes[dataMat[,2] != 0 & dataMat[,3] != 0 & dataMat[,4] != 0] = col = "black"
####################################################################
# The color of the heatmap itself gets the colors; purple, lightgrey and green.
hmcol <- colorRampPalette(c("purple", "lightgrey", "green")) (n=99)
####################################################################
#                         LogFC heatmap                            #
####################################################################
# The real data is obtained with the use of the M3 dataset
# The values of the obtained data will be created with the use of the 
# raw information - the mean of the row. This creates the logFC.
logFCVal <- M3[match(rownamesData, rownames(M3)),1:12] - apply(M3[match(rownamesData, rownames(M3)),1:12], 1, mean)
# The colors for the genes are also added into the matrix
logFCVal <- apply(cbind(logFCVal, col_genes), 2, function(x) unname(unlist(x)))
# The rownames will be the genes.
rownames(logFCVal) <- rownames(M3[match(rownamesData, rownames(M3)),1:12] )
# An ordering is done on the group were the genes come from and the highest 
# expression of the HET6 mice.
logFCVal <- logFCVal[order(as.numeric(as.character(logFCVal[,12])), decreasing = T),]
logFCVal <- logFCVal[order(factor(logFCVal[,13], levels = c("green", "palegreen", "grey", "orchid", "purple", "black"))),]
# The colors for the genes are cut off into another matrix.
logFCResults <- logFCVal[,1:12]
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
