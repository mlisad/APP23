####################################################################
# Author    : M. Dubbelaar
# Date      : 08-okt-2015
# File Name : LFCseqVisualisation.R
# Purpose   : To visualise the APP data with the package LFCseq.
# Used files: LFCseqR.R
#             LFCSeqR_helper.R
#             LoadingAppfile.R
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/Desktop/Onderzoek-APP23_RNASEQ/CreatedData/LfcSeq/")
source("Code/LFCseqR.R")
source("Code/LFCseqR_helper.R")
#install.packages("RColorBrewer")
library("RColorBrewer")
library("gplots")
####################################################################
source("../EdgeR/Code/EdgeRFunctions.R")
M1 <- getData("/Users//mldubbelaar/Downloads/expression_table02.genelevel.GRCm38.v76.htseq.txt.table")
targets <- getTarget("/Users//mldubbelaar/APP23/Targets.csv")
source("../plotColors.R")
####################################################################
#                       Using the main function                    #
####################################################################
# The data returned by the function LFCseq are the p values.
# The data is calculated with the use of three differen normalisation
# methods. First the data is filtered to reduce any background noise.
# The depth will be estimated next with the use of three methods
# and least an estimation will be made for the non DE data of
# each feature.
# The LfcSeq function uses all of the function which are found in the
# source LFCseqR_helper.R
hmcol <- colorRampPalette(brewer.pal(9,"YlGnBu"))(100)

saveFoundData <- function(data, samples, norm, mapName) {
  normData <- LFCseq(data, condsAB = targets$Conditie, norm.method = norm)
  dataMat <- data.frame(row.names = rownames(data), normData)
  unique <- which(dataMat < 0.05)
  unique <- data[unique,]
  colnames(unique) <- samples
  pdf(paste("Plots/", mapName, norm, "Heatmap.pdf", sep=""))
  heatmap.2(as.matrix(unique), cexRow = 0.0001,  col = hmcol, trace = "none", scale = "row", main = paste("Data ", data))
  dev.off()
  write.table(rownames(unique), file = paste("/home/mdubbelaar/Desktop/APP23_results/LfcSeq/Made_Documents/", mapName, norm, "Genes.txt", sep=""), row.names = F, col.names=F, sep=", ", quote = F)
}

normData <- LFCseq(M1, condsAB = targets$Conditie, norm.method = "rpm")
dataMat <- data.frame(row.names = rownames(M1), normData)
unique <- which(dataMat < 0.05)
unique <- M1[unique,]
colnames(unique) <- targets$Samples

mostExpr <- order(rowMeans(unique), decreasing = T)[1:30]
mostExprInfo <- unique[mostExpr,]
dingen <- apply(mostExprInfo, 2, sum)
barplot(dingen, col=c(rep((col=rainbow(265)[20]), 3), rep((col=rainbow(265)[145]), 3)))

saveFoundData(M1, targets$Samples, "rpm", "All/")
saveFoundData(M1, targets$Samples, "deseq","All/")
saveFoundData(M1, targets$Samples, "npseq", "All/")

saveFoundData(M1[,7:24], targets$Samples[7:24], "rpm", "old_Mice/")
saveFoundData(M1[,7:24], targets$Samples[7:24], "deseq","old_Mice/")
saveFoundData(M1[,7:24], targets$Samples[7:24], "npseq", "old_Mice/")

saveFoundData(M1[,1:12], targets$Samples[1:12], "rpm", "young_Mice/")
saveFoundData(M1[,1:12], targets$Samples[1:12], "deseq","young_Mice/")
saveFoundData(M1[,1:12], targets$Samples[1:12], "npseq", "young_Mice/")
