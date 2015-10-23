####################################################################
# Author    : M. Dubbelaar
# Date      : 28-sept-2015
# File Name : LimmaDispersions.R
# Purpose   : Calculates the dispersions of each comparison
####################################################################
# The design is used for the calculation of the estimates.
design <- model.matrix(~0+factor(targets$Conditie), data = dge$samples)
# Voom contains several build-in procedures to estimate the main-variance relationship.
v <- voom(dge, design, plot = T)
# lmfit fits the linear model.
fit <- lmFit(v, design)

# The cont.matrix consists of all comaprisons which need to be made.
# The DE will be defined and the contrast is shown behind it.
cont.matrix <- cbind("6_8W.WTvs6_8.W.HET"=c(1,0,0,0,-1,0,0,0),
                     "12M.WTvs12M.HET"=c(0,1,0,0,0,-1,0,0),
                     "18M.WTvs18M.HET"=c(0,0,1,0,0,0,-1,0),
                     "24M.WTvs24M.HET"=c(0,0,0,1,0,0,0,-1),
                     "12M.HETvs6_8W.HET"=c(1,-1,0,0,0,0,0,0),
                     "18M.HETvs12M.HET"=c(0,1,-1,0,0,0,0,0),
                     
                     "24M.HETvs18M.HET"= c(0,0,1,-1,0,0,0,0),
                     "12M.WTvs6_8W.WT"=c(0,0,0,0,1,-1,0,0),
                     "18M.WTvs12M.WT"=c(0,0,0,0,0,1,-1,0),
                     "24M.WTvs18M.WT"=c(0,0,0,0,0,0,1,-1),
                     "24M.HETvs6_8W.HET"=c(1,0,0,-1,0,0,0,0),
                     "24M.WTvs6_8W.WT"=c(0,0,0,0,1,0,0,-1))
########################################################################################################################################
# The fit is made with the use of the eBayes function (note not ebayes this causes errors).
fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
# Data is stored within the tables for further anaylsis.
# The coef is defined, n = describes the amount of genes which can be saved within the table and the sorting is none
# the adjusted.method is not used, because the default is "BH".
limma.table1 <- topTable(fit2, coef = "6_8W.WTvs6_8.W.HET", n=Inf, sort.by="none")
limma.table2 <- topTable(fit2, coef = "12M.WTvs12M.HET", n=Inf, sort.by="none")
limma.table3 <- topTable(fit2, coef = "18M.WTvs18M.HET", n=Inf, sort.by="none")
limma.table4 <- topTable(fit2, coef = "24M.WTvs24M.HET", n=Inf, sort.by="none")
limma.table5 <- topTable(fit2, coef = "12M.HETvs6_8W.HET", n=Inf, sort.by="none")
limma.table6 <- topTable(fit2, coef = "18M.HETvs12M.HET", n=Inf, sort.by="none")
limma.table7 <- topTable(fit2, coef = "24M.HETvs18M.HET", n=Inf, sort.by="none")
limma.table8 <- topTable(fit2, coef = "12M.WTvs6_8W.WT", n=Inf, sort.by="none")
limma.table9 <- topTable(fit2, coef = "18M.WTvs12M.WT", n=Inf, sort.by="none")
limma.table10 <- topTable(fit2, coef = "24M.WTvs18M.WT", n=Inf, sort.by="none")
limma.table11 <- topTable(fit2, coef = "24M.HETvs6_8W.HET", n=Inf, sort.by="none")
limma.table12 <- topTable(fit2, coef = "24M.WTvs6_8W.WT", n=Inf, sort.by="none")
########################################################################################################################################
# The results from the toptables are checked with the use of the own-made filtergenes function
# This function returns the data which contain a p-value below the 0.05
toptable1.results <- filterGenesWithLimma(limma.table1)
toptable2.results <- filterGenesWithLimma(limma.table2)
toptable3.results <- filterGenesWithLimma(limma.table3)
toptable4.results <- filterGenesWithLimma(limma.table4)
toptable5.results <- filterGenesWithLimma(limma.table5)
toptable6.results <- filterGenesWithLimma(limma.table6)
toptable7.results <- filterGenesWithLimma(limma.table7)
toptable8.results <- filterGenesWithLimma(limma.table8)
toptable9.results <- filterGenesWithLimma(limma.table9)
toptable10.results <- filterGenesWithLimma(limma.table10)
toptable11.results <- filterGenesWithLimma(limma.table11)
toptable12.results <- filterGenesWithLimma(limma.table12)