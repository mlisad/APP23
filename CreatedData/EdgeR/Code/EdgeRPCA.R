####################################################################
# Author    : M. Dubbelaar
# Date      : 15-sept-2015
# File Name : Creating_PCA.r
# Purpose   : Creates a PCA 3D plot (saved in 2D)
####################################################################
#              Installing all of the necessary packages            #
####################################################################
#install.packages("rgl")
#install.packages("gplots")
library(rgl)

pcaPlot <- function(pathway) {
####################################################################
#                       Preparing data for PCA                     #
####################################################################
# The dataset is rotated and after this the coloms will be filtered
# When the sum of the data is a negative number, this colomn will be deleted.
tset <- t(M2)
tset <- tset[,apply(tset, 2, sum) > 0]
pca <- prcomp(tset, cor=T, scores= T, scale=T) 
####################################################################
#                           Creating PCA                           #
####################################################################
# open3d is used to save the data
open3d()
# The PCA plot is made with the use of the amount of samples.
plot3d(pca$x, col=col_cell_age, size = "5")
# rgl.postscript writes the current figure into a pdf (so you can rotate the plot before saving). 
rgl.postscript(paste(pathway, "Plots/PCA.pdf", sep=""), "pdf")
}