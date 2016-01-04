####################################################################
# Author    : M. Dubbelaar
# Date      : 28-sept-2015
# File Name : LimmaLinearTime.R
# Purpose   : Calculates with the use of multiple designs the linear
#             genotype agains age, the main effect of age and the 
#             main effect of genotype
####################################################################
#                            All results                           #
####################################################################
# The group will be determined and must be releved to use the
# WT as the intercept. The data can be used by the function
# with the use of the group and the design
genotypeGroup <- factor(targets$Genotype)
genotypeGroup <- relevel(genotypeGroup, "WT")
design <- model.matrix(~genotypeGroup*targets$Age)

calculateEffectsLimma(dge, design, paste(resultPathway, "Made_Documents/All_ages/main_genotype_result.txt", sep=""), paste(resultPathway,"Made_Documents/All_ages/main_age_result.txt", sep="")
                        , paste(resultPathway, "Made_Documents/All_ages/interaction_result.txt", sep=""), paste(resultPathway, "Plots/linear_time_heatmaps/All_ages_heatmaps.pdf", sep=""))
####################################################################
#                   All results without 6-8 weeks                  #
####################################################################
# The group of the old mice doesn't contain the mice with an 
# age of 6-8 weeks (2 months). The releveling takes place to make
# sure that the WT is used as the intercept.
genotypeOldMiceGroup <- factor(targets$Genotype[targets$Age != 2])
genotypeOldMiceGroup <- relevel(genotypeOldMiceGroup, "WT")
oldMiceDesign <- model.matrix(~genotypeOldMiceGroup*targets$Age[targets$Age != 2])
calculateEffectsLimma(dge[,7:24], oldMiceDesign, paste(resultPathway ,"Made_Documents/6-18-24M_old_mice/main_genotype_result.txt", sep=""), paste(resultPathway, "Made_Documents/6-18-24M_old_mice/main_age_result.txt", sep="")
                        , paste(resultPathway, "Made_Documents/6-18-24M_old_mice/interaction_result.txt", sep=""), paste(resultPathway, "Plots/linear_time_heatmaps/old_mice_heatmaps.pdf", sep=""))
####################################################################
#             All results with 6-8 weeks and 12 months             #
####################################################################
# The group of the young mice consist of the 2 and 12 months old
# mice only. The factor is filled in by hand to make sure that the
# data doesn't create any errors. Even here the group must be
# releveled to make sure that the WT will be used as the intercept.
genotypeYoungMiceGroup <- factor(rep(c(rep("WT",3), rep("HET",3)), 2))
genotypeYoungMiceGroup <- relevel(genotypeYoungMiceGroup, "WT")
youngMiceDesign <- model.matrix(~genotypeYoungMiceGroup*c(targets$Age[targets$Age == "2" ], targets$Age[targets$Age == "6"]))
calculateEffectsLimma(dge[,1:12], youngMiceDesign, paste(resultPathway, "Made_Documents/2M-6M_old_mice/main_genotype_result.txt", sep=""), paste(resultPathway, "Made_Documents/2M-6M_old_mice/main_age_result.txt", sep="")
                        , paste(resultPathway, "Made_Documents/2M-6M_old_mice/interaction_result.txt", sep=""), paste(resultPathway, "Plots/linear_time_heatmaps/young_mice_heatmaps.pdf", sep=""))