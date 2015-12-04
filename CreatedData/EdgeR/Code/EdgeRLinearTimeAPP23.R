####################################################################
# Author    : M. Dubbelaar
# Date      : 23-sept-2015
# File Name : calculate_linear_timepoints.R
# Purpose   : Calculates with the use of multiple designs the linear
#             genotype agains age, the main effect of age and the 
#             main effect of genotype
####################################################################
#               Results with a design Genotype * Time              #
####################################################################
# The Group is used for the design and the releveling is necessary 
# to use the WT as the intercept.
# All of the data will be compared to the Wild Type instead of the APP mice.
Group <- factor(targets$Genotype)
Group <- relevel(Group, "WT")
# The timpoints are the different ages (2, 6, 18 and 24 months)
design.timepoints <- model.matrix(~Group*targets$Age)
dge.timepoints <- estimateGLMCommonDisp(dge, design.timepoints)
dge.timepoints <- estimateGLMTrendedDisp(dge.timepoints, design.timepoints)
calculateMaineffectsInteraction(dge.timepoints, design.timepoints, "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Made_Documents/All_ages/", "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/linear_time_heatmaps/all_ages_heatmaps.pdf", c(1:24))

####################################################################
#    Results with a design Genotype * Time without 6-8 weeks       #
####################################################################
# The mice with an age of 2 months are excluded from the next design.
# A releveling is used again to use the WT data as an intercept.
olderGroup <- factor(targets$Genotype[targets$Age != 2])
olderGroup <- relevel(olderGroup, "WT")
# The matrix is made without the 2 months old mice
designOlderMice <- model.matrix(~olderGroup*targets$Age[targets$Age != "2"])
dgeOlderMice <- estimateGLMCommonDisp(dge[,7:24], designOlderMice)
dgeOlderMice <- estimateGLMTrendedDisp(dgeOlderMice, designOlderMice)
resultsOlderMice <- calculateMaineffectsInteraction(dgeOlderMice, designOlderMice, "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Made_Documents/6-18-24M_old_mice/", "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/linear_time_heatmaps/old_mice_heatmaps.pdf", c(7:24))

saveInfoDE("APP23", resultsOlderMice, "DifferentialGenesMainGenotypeOldMice.txt", "DifferentialGenesMainAgeOldMice.txt", "DifferentialGenesLinearOldMice.txt")
plotMostExpr(resultsOlderMice, which(resultsOlderMice[[2]][[1]]$FDR < 0.01 & resultsOlderMice[[2]][[1]]$logFC > .075), which(resultsOlderMice[[3]][[1]]$FDR < 0.01 & resultsOlderMice[[3]][[1]]$logFC > .05),
             7:24, "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/AgeAndLinearOld.pdf")
####################################################################
#   Results with a design Genotype * Time without 18 & 24 months   #
####################################################################
# The rep of the data is made and is releveled to make the WT the intercept.
youngerGroup <- factor(rep(c(rep("WT",3), rep("HET",3)),2))
youngerGroup <- relevel(youngerGroup, "WT")
# Only the data of the mice with an age of 2 and 12 (the youngest mice) are used.
designYoungerMice <- model.matrix(~youngerGroup*c(targets$Age[targets$Age == "2" ], targets$Age[targets$Age == "6"]))
dgeYoungerMice <- estimateGLMCommonDisp(dge[,1:12], designYoungerMice)
dgeYoungerMice <- estimateGLMTrendedDisp(dgeYoungerMice, designYoungerMice)
resultsYoungerMice <- calculateMaineffectsInteraction(dgeYoungerMice, designYoungerMice, "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Made_Documents/2M-6M_old_mice/","/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/linear_time_heatmaps/young_mice_heatmaps.pdf", c(1:12))

saveInfoDE("APP23", resultsYoungerMice,"DifferentialGenesMainGenotypeYoungMice.txt", "DifferentialGenesMainAgeYoungMice.txt", "DifferentialGenesLinearYoungMice.txt")
plotMostExpr(resultsYoungerMice, which(resultsYoungerMice[[2]][[1]]$FDR < 0.01 & resultsYoungerMice[[2]][[1]]$logFC > .15), which(resultsYoungerMice[[3]][[1]]$FDR < 0.05),
             1:12, "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/AgeAndLinearYoung.pdf")
