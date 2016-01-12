####################################################################
# Author    : M. Dubbelaar
# Date      : 16-nov-2015
# File Name : EdgeRScatterplots.R
# Purpose   : Created scatterplots/vulcano plots for the different
#             comparisons of the APP23 dataset.
####################################################################
#                       Create vulcanoplot data                    #
####################################################################
# These results contain the different genes that are calculated with 
# the different pvals and logfc values.
results1 <- calFDRandLogFC(toptable1[[1]], 0.05, 0.05, 0.01, 0.2, 0.3)
results2 <- calFDRandLogFC(toptable2[[1]], 0.05, 0.05, 0.01, 0.2, 0.3)
results3 <- calFDRandLogFC(toptable3[[1]], 0.05, 0.01, 0.005, 0.5, 1.5)
results4 <- calFDRandLogFC(toptable4[[1]], 0.05, 0.01, 0.005, 0.7, 1.7)
results5 <- calFDRandLogFC(toptable5[[1]], 0.05, 0.01, 0.005, 0.8, 1.8)
results6 <- calFDRandLogFC(toptable6[[1]], 0.05, 0.01, 0.005, 1, 2)
results7 <- calFDRandLogFC(toptable7[[1]], 0.05, 0.01, 0.005, 1, 2)
results8 <- calFDRandLogFC(toptable8[[1]], 0.05, 0.01, 0.005, 1, 2)
results9 <- calFDRandLogFC(toptable9[[1]], 0.05, 0.01, 0.005, 0.4, 0.8)
results10 <- calFDRandLogFC(toptable10[[1]], 0.05, 0.01, 0.005, 0.8, 1.8)
results11 <- calFDRandLogFC(toptable11[[1]], 0.05, 0.01, 0.005, 1.5, 3.1)
results12 <- calFDRandLogFC(toptable12[[1]], 0.05, 0.01, 0.005, 1.5, 3)
####################################################################
#                       Create images from data                    #
####################################################################
# The code below make different files with the different comparisons 
# of the APP23 dataset.
# The first set creates vulcanoplots with the DEG of the genotypes.
pdf(paste(resultPathway, "Plots/Vulcano_plot_WT-HET.pdf", sep=""), height = 10, width= 10)
par(mfrow=c(4,4))
makeVulcanoPlot(toptable1, results1, "6-8 weeks WT vs 6-8 weeks HET")
makeVulcanoPlot(toptable2, results2, "6 months WT vs 6 months HET")
makeVulcanoPlot(toptable3, results3, "18 months WT vs 18 months HET")
makeVulcanoPlot(toptable4, results4, "24 months WT vs 24 months HET")
dev.off()
# The second set creates vulcanoplots with DEG among the different ages in the HET mice.
pdf(paste(resultPathway, "Plots/Vulcano_plot_HET.pdf", sep=""), height = 10, width= 10)
par(mfrow=c(1,3))
makeVulcanoPlot(toptable5, results5, "6 months HET vs 6-8 weeks HET")
makeVulcanoPlot(toptable6, results6, "18 months HET vs 6 months HET")
makeVulcanoPlot(toptable7, results7, "24 months HET vs 18 months HET")
dev.off()
# The third set creates vulcanoplots with DEG among the different ages in the WT mice.
pdf(paste(resultPathway, "Plots/Vulcano_plot_WT.pdf", sep=""), height = 10, width= 10)
par(mfrow=c(1,3))
makeVulcanoPlot(toptable8, results8, "6 months WT vs 6-8 weeks WT")
makeVulcanoPlot(toptable9, results9, "18 months WT vs 6 months WT")
makeVulcanoPlot(toptable10, results10, "24 months WT vs 18 months WT")
dev.off()
# The fourth set creates the vulcanoplots with DEG between the oldest and the youngest mice
# for both the genotypes.
pdf(paste(resultPathway, "Plots/Vulcano_plot_old-young.pdf", sep=""), height = 10, width= 10)
par(mfrow=c(1,2))
makeVulcanoPlot(toptable11, results11, "24 months HET vs 6-8 weeks HET")
makeVulcanoPlot(toptable11, results12, "24 months WT vs 6-8 weeks WT")
par(mfrow=c(1,1))
dev.off()
