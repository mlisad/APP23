####################################################################
# Author    : M. Dubbelaar
# Date      : 23-sept-2015
# File Name : calculate_linear_timepoints.R
# Purpose   : Calculates with the use of multiple designs the linear
#             genotype agains age, the main effect of age and the 
#             main effect of genotype
####################################################################
#                Function for calculating the effects              #
####################################################################
calculateMaineffectsInteraction <- function(dge, design, pathwayDoc, pathwayPlot, length) {
  # This function makes sure that the main effect of the ages, the main effect of the genotype and the interaction model
  # of each dataset is calculated.
  # The dge and the design are necessary items to calculate the these effects.
  # glmFit conducts statistical tests to fit the negative binomial generalized linear model.
  # The function glmLRT, tests the likelihood of the coefficients.
  # topTags orders the data by the ranking of the p value ord the logFC.
  # The lasts steps contains the filtering of the genes with a FDR < 0.05, the visualisation of these genes and 
  # saving these unique gene names into a file.

  fit <- glmFit(dge, design)
  # With the use of different coefs, different genewise statistical tests are made.
  # The first one compares the WT against APP (main effect of genotype)
  lrt1 <- glmLRT(fit, coef=2)
  # The second comparison checks the main effect of age
  lrt2 <- glmLRT(fit, coef=3)
  # The third and last comparison checks the linear interaction.
  lrt3 <- glmLRT(fit, coef=4)
  
  # The 3 different statistical tests are extracted in a dataframe
  # These results will be checked within the own-made function filtergenes to check if they meet the requirements.
  Toptable1 <- topTags(lrt1, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  Toptable2 <- topTags(lrt2, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  Toptable3 <- topTags(lrt3, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  
  Toptable1.results <- filterGenes(Toptable1)
  geneInfo1 <- Toptable1.results[[2]]
  Toptable1.results <- Toptable1.results[[1]]

  Toptable2.results <- filterGenes(Toptable2)
  geneInfo2 <- Toptable2.results[[2]]
  Toptable2.results <- Toptable2.results[[1]]
  
  Toptable3.results <- filterGenes(Toptable3)
  geneInfo3 <- Toptable3.results[[2]]
  Toptable3.results <- Toptable3.results[[1]]
  
  # The last step is to plot the genes with their information in a heatmap.
  pdf(pathwayPlot) 
  heatmap.2(M2[match(rownames(Toptable1.results), rownames(M2)), length], ColSideColors= col_cell_age[length], cexRow = 0.01, trace = "none", scale = "row", main="Main effect genotype")
  heatmap.2(M2[match(rownames(Toptable2.results), rownames(M2)), length], ColSideColors = col_cell_age[length], cexRow = 0.01, trace = "none", scale = "row", main="Main effect age")
  heatmap.2(M2[match(rownames(Toptable3.results), rownames(M2)), length], ColSideColors = col_cell_age[length], cexRow = 0.01, trace = "none", scale = "row", main="Interaction effect")
  dev.off()
  
  # The unique genes are saved within a table so these gene names can be used for futher analysis. 
  write.table(rownames(Toptable1.results), paste(pathwayDoc, "main_genotype_result.txt", sep = ""), row.names = F, col.names=F, eol=",\n", quote = F)
  write.table(rownames(Toptable2.results), paste(pathwayDoc, "main_age_result.txt", sep = ""), row.names = F, col.names=F, eol=",\n", quote = F)
  write.table(rownames(Toptable3.results), paste(pathwayDoc, "interaction_result.txt", sep = ""), row.names = F, col.names=F, eol=",\n", quote = F)

  DE_Expression <- cbind(rownames(Toptable1[[1]]), Toptable1[[1]]$logFC, Toptable1[[1]]$FDR, Toptable2[[1]]$logFC, Toptable2[[1]]$FDR, Toptable3[[1]]$logFC, Toptable3[[1]]$FDR, BioM[,3:4])
  write.table(DE_Expression , paste(pathwayDoc, "LinearTimeDE.txt", sep=""), row.names = F,  col.names = c("Genes", "Main-effect Genotype (logFC)", 
                            "Main-effect Gentyope (FDR)", "Main-effect Age (logFC)", "Main-effect Age (FDR)", "Linear Effect (logFC)", 
                            "Linear Effect (FDR)", "Gene Symbol", "Gene Description"), sep = "\t")
  data <- list(Toptable2.results, Toptable3.results, geneInfo2, geneInfo3)
  }
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

saveInfoDE(resultsOlderMice, "DifferentialGenesMainAgeOldMice.txt", "DifferentialGenesLinearOldMice.txt")
plotMostExpr(resultsOlderMice, which(resultsOlderMice[[1]][[1]]$FDR < 0.01 & resultsOlderMice[[1]][[1]]$logFC > .075), which(resultsOlderMice[[2]][[1]]$FDR < 0.01 & resultsOlderMice[[2]][[1]]$logFC > .05),
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

saveInfoDE(resultsOlderMice, "DifferentialGenesMainAgeYoungMice.txt", "DifferentialGenesLinearYoungMice.txt")
plotMostExpr(resultsYoungerMice, which(resultsYoungerMice[[1]][[1]]$FDR < 0.01 & resultsYoungerMice[[1]][[1]]$logFC > .15), which(resultsYoungerMice[[2]][[1]]$FDR < 0.05),
             1:12, "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/30mostexpressedGenes/AgeAndLinearYoung.pdf")
