####################################################################
# Author    : M. Dubbelaar
# Date      : 08-okt-2015
# File Name : LFCseqVisualisation.R
# Purpose   : To visualise the APP data with the package LFCseq.
# Used files: LFCseqR.R
#             LFCSeqR_helper.R
#             EdgeRFunctions.R
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/APP23/CreatedData/LfcSeq/")
#install.packages("RColorBrewer")
#install.packages("gplots")
library("RColorBrewer")
library("gplots")
source("Code/LFCseqR.R")
source("Code/LFCseqR_helper.R")
resultPathway <- "/home/mdubbelaar/Desktop/APP23_results/LfcSeq/"
####################################################################
#                            Loading data                          #
####################################################################
source("../EdgeR/Code/EdgeRFunctions.R")
M1 <- getData("/media/mdubbelaar/6CEC0BDEEC0BA186/1507_Holtman_RNAseq/run01/results/expression/expressionTable/expression_table02.genelevel.GRCm38.v76.htseq.txt.table")
targets <- getTarget("../../Targets.csv")
####################################################################
#                        Necessary Functions                       #
####################################################################
saveFoundData <- function(data, samples, norm, mapName) {
  # The function LFCseq uses the different functions that are already
  # available in the other files of lfcSeq. This function normalizes the data.
  # The first step is to filter the data to reduce the background noise.
  # The dept is estimated by the method that is chosen, with this dept
  # an estimation of the non DE data is made for each feature.
  # The outcome of the data are p values.
  normData <- LFCseq(data, condsAB = targets$Conditie, norm.method = norm)
  # The ensemble id are put together with the pvalues to create a 
  # data matrix.
  dataMat <- data.frame(row.names = rownames(data), normData)
  # The values below a p value of 0.05 are saved from the M1 matrix
  # and will be used further.
  unique <- which(dataMat < 0.05)
  unique <- data[unique,]
  # The column names will be changes into the sample id's
  colnames(unique) <- samples
  # The data will be saved as a heatmap and the genes found will be saved
  # into a text file.
  pdf(paste("Plots/", mapName, norm, "Heatmap.pdf", sep=""))
  heatmap.2(as.matrix(unique), cexRow = 0.0001,  col = hmcol, trace = "none", scale = "row", main = paste("Data ", data))
  dev.off()
  write.table(rownames(unique), file = paste(resultPathway, "Made_Documents/", mapName, norm, "Genes.txt", sep=""), row.names = F, col.names=F, sep=", ", quote = F)
}
####################################################################
#                       Using the main function                    #
####################################################################
# The vector hmcol contains a colorpattern (yellowGreenBlue) that can 
# be used for the heatmaps.
hmcol <- colorRampPalette(brewer.pal(9,"YlGnBu"))(100)

# Checking all of the different normalization methods for all ages.
saveFoundData(M1, targets$Samples, "rpm", "All/")
saveFoundData(M1, targets$Samples, "deseq","All/")
saveFoundData(M1, targets$Samples, "npseq", "All/")
# Checking all of the different normalization methods for all the 
# old mice (6, 18 and 24 months old mice).
saveFoundData(M1[,7:24], targets$Samples[7:24], "rpm", "old_Mice/")
saveFoundData(M1[,7:24], targets$Samples[7:24], "deseq","old_Mice/")
saveFoundData(M1[,7:24], targets$Samples[7:24], "npseq", "old_Mice/")
# Checking all of the different normalization methods for all the 
# young mice (2 and 6 months old mice).
saveFoundData(M1[,1:12], targets$Samples[1:12], "rpm", "young_Mice/")
saveFoundData(M1[,1:12], targets$Samples[1:12], "deseq","young_Mice/")
saveFoundData(M1[,1:12], targets$Samples[1:12], "npseq", "young_Mice/")
