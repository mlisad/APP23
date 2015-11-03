####################################################################
# Author    : M. Dubbelaar
# Date      : 10-sept-2015
# File Name : loadingAppFile.R
# Purpose   : Loads the data from the APP23 research and the information of the samples.
####################################################################
# The APP23_data contains the probe names and all of the measurement for each sample
APP23_data <- read.delim("/media/mdubbelaar/6CEC0BDEEC0BA186/1507_Holtman_RNAseq/run01/results/expression/expressionTable/expression_table02.genelevel.GRCm38.v76.htseq.txt.table")
APP23_data <- APP23_data[order(colnames(APP23_data), decreasing = F )]

# The information of the APP23_data data is stored within a different file.
# This file contains information like: the sample_nr, the condition, the number of the plate and so on.
targets <- read.table("/home/mdubbelaar/APP23/Targets.csv", sep=":", header = T)
targets <- targets[order(targets$Sample),]
####################################################################

