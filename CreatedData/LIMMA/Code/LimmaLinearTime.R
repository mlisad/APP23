####################################################################
# Author    : M. Dubbelaar
# Date      : 28-sept-2015
# File Name : linear_timepoints_limma.r
# Purpose   : Calculates with the use of multiple designs the linear
#             genotype agains age, the main effect of age and the 
#             main effect of genotype
####################################################################
#                 Function calculate effects limma                 #
####################################################################
calculateEffectsLimma <- function(dge, design, namefile1, namefile2, namefile3, pathway) {
  # This function calculates the main effect of genotype, the main effect of age
  # and the interaction between age and genotype. The function voom transforms the 
  # count data to logCPM (log2-counts per million), estimates the main-variance 
  # relationship. The data can be used for linear modelling after these steps.
  # The fit (calculated with the use of the eBayes function) can be used to
  # tests the likelihood of the coefficients. The genes will be filtered 
  # only genes with an adjusted value below the 0.05 will be saved into the
  # result vector. The saved genes will be used in the heatmap for visualisation
  # and the gene names will be saved within a txt file.
  
  # Voom contains several build-in procedures to estimate the main-variance relationship.
  v <- voom(dge, design)
  # lmfit fits the linear model.
  fit <- lmFit(v, design)
  # The fit is made with the use of the eBayes function (note not ebayes this causes errors).
  fit2 <- eBayes(lmFit(v, design))
  
  # The main effects and the interaction effect are defined.
  mainGenotype <- topTable(fit2, coef=2, n=Inf, sort.by="none")
  mainTime <- topTable(fit2, coef = 3, n=Inf, sort.by="none")
  interactionGenotypeTime <- topTable(fit2, coef = 4,  n=Inf, sort.by="none")
  # and the found genes are filtered for an adjusted p value below the 0.05
  mainGenotype.result <- filterGenesWithLimma(mainGenotype)
  mainTime.result <- filterGenesWithLimma(mainTime)
  interactionGenotypeTime.result <- filterGenesWithLimma(interactionGenotypeTime)
  
  # The last step is to plot the genes with their information in a heatmap.
  # The if - else is necessary because of a lack op interaction genes within one
  # of the age groups.
  if (length(rownames(interactionGenotypeTime.result)) == 0) {
    pdf(pathway)
    heatmap.2(M2[rownames(mainTime.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect age")
    heatmap.2(M2[rownames(mainGenotype.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect genotype")
    dev.off()  
    } else {
    pdf(pathway)
    heatmap.2(M2[rownames(mainTime.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect age")
    heatmap.2(M2[rownames(mainGenotype.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect genotype")
    heatmap.2(M2[rownames(interactionGenotypeTime.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Interaction effect")
    dev.off()
    }
  # The unique genes are saved within a table so these gene names can be used for futher analysis.
  write.table(rownames(mainGenotype.result), namefile1, row.names = F, col.names=F, sep=", ", quote = F)
  write.table(rownames(mainTime.result), namefile2, row.names = F, col.names=F, sep=", ", quote = F)
  write.table(rownames(interactionGenotypeTime.result), namefile3, row.names = F, col.names=F, sep=", ", quote = F)
}
####################################################################
#                            All results                           #
####################################################################
# The group will be determined and must be releved to use the
# WT as the intercept. The data can be used by the function
# with the use of the group and the design
genotypeGroup <- factor(targets$Genotype)
genotypeGroup <- relevel(genotypeGroup, "WT")
design <- model.matrix(~genotypeGroup*targets$Age)

calculateEffectsLimma(dge, design, "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/All_ages/main_genotype_result.txt", "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/All_ages/main_age_result.txt"
                        , "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/All_ages/interaction_result.txt", "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Plots/linear_time_heatmaps/All_ages_heatmaps.pdf" )
####################################################################
#                   All results without 6-8 weeks                  #
####################################################################
# The group of the old mice doesn't contain the mice with an 
# age of 6-8 weeks (2 months). The releveling takes place to make
# sure that the WT is used as the intercept.
genotypeOldMiceGroup <- factor(targets$Genotype[targets$Age != 2])
genotypeOldMiceGroup <- relevel(genotypeOldMiceGroup, "WT")
oldMiceDesign <- model.matrix(~genotypeOldMiceGroup*targets$Age[targets$Age != 2])
calculateEffectsLimma(dge[,7:24], oldMiceDesign, "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/6-18-24M_old_mice/main_genotype_result.txt", "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/6-18-24M_old_mice/main_age_result.txt"
                        , "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/6-18-24M_old_mice/interaction_result.txt", "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Plots/linear_time_heatmaps/old_mice_heatmaps.pdf" )
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
calculateEffectsLimma(dge[,1:12], youngMiceDesign, "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/2M-6M_old_mice/main_genotype_result.txt", "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/2M-6M_old_mice/main_age_result.txt"
                        , "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Made_Documents/2M-6M_old_mice/interaction_result.txt", "/home/mdubbelaar/Desktop/APP23_results/LIMMA/Plots/linear_time_heatmaps/young_mice_heatmaps.pdf" )
