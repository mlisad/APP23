####################################################################
# Author    : M. Dubbelaar
# Date      : 03-dec-2015
# File Name : EdgeRLinearTimeCKp25.R
# Purpose   : Calculates with the use of a linear design the 
#             genotype agains age, the main effect of age and the 
#             main effect of genotype.
####################################################################
#               Results with a design Genotype * Time              #
####################################################################
# The Group is used for the design and the releveling is necessary 
# to use the WT as the intercept.
# All of the data will be compared to the Wild Type instead of the Ckp25 mice.
Group <- factor(targets$Genotype)
Group <- relevel(Group, "WT")
# The timpoints are the different ages (2 and 6 weeks)
design.timepoints <- model.matrix(~Group*targets$Age)
dge.timepoints <- estimateGLMCommonDisp(dge, design.timepoints)
dge.timepoints <- estimateGLMTrendedDisp(dge.timepoints, design.timepoints)
resultsCkp25 <- calculateMaineffectsInteraction(dge.timepoints, design.timepoints, "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/", "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/GLM_heatmaps.pdf", c(1:12))
head(resultsCkp25[[1]][[1]])
saveInfoDE("CKp25", resultsCkp25, "DifferentialGenesMainGenotype.txt", "DifferentialGenesMainAge.txt", "DifferentialGenesLinear.txt")


