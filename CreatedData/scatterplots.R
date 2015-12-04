
calFDRandLogFC <- function (tablename, pval1, pval2, pval3, logfc, logfc2) {
  Toptable1 <- tablename[abs(as.data.frame(tablename[,5])) < pval1,] 
  Toptable2 <- tablename[abs(as.data.frame(tablename[,5])) < pval2 & abs(as.data.frame(tablename[,1])) > logfc,]
  Toptable3 <- tablename[abs(as.data.frame(tablename[,5])) < pval3 & abs(as.data.frame(tablename[,1])) > logfc2,]
  Genes1 <- BioM[abs(as.data.frame(tablename[,5])) < pval3 &  abs(as.data.frame(tablename[,1])) > logfc2, 3]
  results <- list(Toptable1, Toptable2, Toptable3, Genes1)
  return (results)
}

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
############################################################
changeFDR <- function(table) {
  pval <- table[[1]][1]
  pval[pval < 1E-70,] <- 1E-70 
  return (pval)
}

makeplot <- function(table, result, name) {
  plot(as.data.frame(table[[1]][,1])[[1]], -log(as.data.frame(changeFDR(table[[1]][,5]))[[1]], 10), main=name, pch=20, xlab="log2FC)", ylab="-log10 (FDR p-value)")
  points(as.data.frame(result[[1]][,1])[[1]], -log(as.data.frame(changeFDR(result[[1]][,5]))[[1]], 10), cex = 1, col= "darkred", pch= 20)
  points(as.data.frame(result[[2]][,1])[[1]], -log(as.data.frame(changeFDR(result[[2]][,5]))[[1]], 10), cex = 1, col= "red", pch= 20)
  text(as.data.frame(result[[3]][,1])[[1]], -log(as.data.frame(changeFDR(result[[3]][,5]))[[1]], 10), result[[4]], cex = 0.75, col= "black")
}
############################################################
pdf("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/Vulcano_plot_WT-HET.pdf", height = 10, width= 10) 
par(mfrow=c(4,4))
makeplot(toptable1, results1, "6-8 weeks WT vs 6-8 weeks HET")
makeplot(toptable2, results2, "6 months WT vs 6 months HET")
makeplot(toptable3, results3, "18 months WT vs 18 months HET")
makeplot(toptable4, results4, "24 months WT vs 24 months HET")
dev.off()

pdf("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/Vulcano_plot_HET.pdf", height = 10, width= 10) 
par(mfrow=c(1,3))
makeplot(toptable5, results5, "6 months HET vs 6-8 weeks HET")
makeplot(toptable6, results6, "18 months HET vs 6 months HET")
makeplot(toptable7, results7, "24 months HET vs 18 months HET")
dev.off()

pdf("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/Vulcano_plot_WT.pdf", height = 10, width= 10) 
par(mfrow=c(1,3))
makeplot(toptable8, results8, "6 months WT vs 6-8 weeks WT")
makeplot(toptable9, results9, "18 months WT vs 6 months WT")
makeplot(toptable10, results10, "24 months WT vs 18 months WT")
dev.off()

pdf("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/Vulcano_plot_old-young.pdf", height = 10, width= 10) 
par(mfrow=c(1,2))
makeplot(toptable11, results11, "24 months HET vs 6-8 weeks HET")
makeplot(toptable11, results12, "24 months WT vs 6-8 weeks WT")
par(mfrow=c(1,1))
dev.off()

pdf("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/Vulcano_plot_24Mwt-24Mhet.pdf", height = 10, width= 10) 
makeplot(toptable4, results4, "24 months WT vs 24 months HET")
dev.off()
