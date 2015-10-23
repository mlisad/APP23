####################################################################
# Author    : M. Dubbelaar
# Date      : 14-sept-2015
# File Name : estimating_dispersions.r
# Purpose   : 
# Used Files: Visualization_APP23_data.R
####################################################################
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)

#my.contrast <- makeContrasts(
#  All.minus6_8weeks = (-WT.24-WT.18-WT.12)/3-(-HET.24-HET.18-HET.12)/3,
#  WT.HET12_18 = (-WT.18-WT.12)/2-(-HET.18-HET.12)/2,
#  WT.HET18_24 = (-WT.24-WT.18)/2-(-HET.24-HET.18)/2,
#  WT18_24 = WT.24-WT.18,
#  HET18_24 = HET.24-HET.18,
#  All.minus12months = (-WT.24-WT.18-WT.06_08)/3-(-HET.24-HET.18-HET.06_08)/3,
#  WT.HET6.8_18 = (-WT.18-WT.06_08)/2-(-HET.18-HET.06_08)/2,
#  WT.HET6.8_24 = (-WT.24-WT.06_08)/2-(-HET.24-HET.06_08)/2,
#  All =(-WT.24-WT.18-WT.12-WT.06_08)/4-(-HET.24-HET.18-HET.12-HET.06_08)/4,
#  levels=design
#)

# Vergelijkt alles met elkaar
lrt01 <- glmLRT(fit, contrast = c(0.25,0.25,0.25,0.25,-0.25,-0.25,-0.25,-0.25))
#lrt1 <- glmLRT(fit, contrast = my.contrast[, "All"])
# Vergelijkt alles behalve de 6/8 weken
lrt02 <- glmLRT(fit, contrast = c(0,0.33,0.33,0.34,0,-0.33,-0.33,-0.34))
# Vergelijkt de 12 en 18 maanden
lrt03 <- glmLRT(fit, contrast = c(0,0.5,0.5,0,0,-0.5,-0.5,0))
# Vergelijkt de 18 en de 24 maanden
lrt04 <- glmLRT(fit, contrast = c(0,0,0.5,0.5,0,0,-0.5,-0.5))
#wt18-24
lrt05 <- glmLRT(fit, contrast = c(0,0,0,0,0,0,-1,1))
#het18-24
lrt06 <- glmLRT(fit, contrast = c(0,0,-1,1,0,0,0,0))
# Vergelijkt alles behalve de 12 maanden
lrt07 <- glmLRT(fit, contrast = c(0.34,0,0.33,0.33,-0.34,0,-0.33,-0.33))
# Vergelijkt 6-8 weken en 18 maand
lrt08 <- glmLRT(fit, contrast = c(0.5,0,0.5,0,-0.5,0,0,-0.5))
# Vergelijkt 6-8 weken en 24 maand
lrt09 <- glmLRT(fit, contrast = c(0.5,0,0,0.5,-0.5,0,0,-0.5))

#het18-wt18
#lrt07 <- glmLRT(fit, contrast = c(0,0,1,0,0,0,-1,0))
#het24-wt24
#lrt08 <- glmLRT(fit, contrast = c(0,0,0,1,0,0,0,-1))
#het18-wt24
#lrt09 <- glmLRT(fit, contrast = c(0,0,1,0,0,0,0,-1))

# Probeer deze lrt's te verminderen
#lrt01 <- glmLRT(fit, contrast = c(-1,0,0,0,1,0,0,0))
#lrt02 <- glmLRT(fit, contrast = c(0,-1,0,0,0,1,0,0))
#lrt03 <- glmLRT(fit, contrast = c(0,0,-1,0,0,0,1,0))
#lrt04 <- glmLRT(fit, contrast = c(0,0,0,-1,0,0,0,1))
#lrt05 <- glmLRT(fit, contrast = c(-1,1,0,0,0,0,0,0))
#lrt06 <- glmLRT(fit, contrast = c(-1,0,1,0,0,0,0,0))
#lrt07 <- glmLRT(fit, contrast = c(-1,0,0,1,0,0,0,0))
#lrt08 <- glmLRT(fit, contrast = c(0,-1,1,0,0,0,0,0))
#lrt09 <- glmLRT(fit, contrast = c(0,-1,0,1,0,0,0,0))
#lrt10 <- glmLRT(fit, contrast = c(0,0,-1,1,0,0,0,0))
#lrt11 <- glmLRT(fit, contrast = c(0,0,0,0,-1,1,0,0))
#lrt12 <- glmLRT(fit, contrast = c(0,0,0,0,-1,0,1,0))
#lrt13 <- glmLRT(fit, contrast = c(0,0,0,0,-1,0,0,1))
#lrt14 <- glmLRT(fit, contrast = c(0,0,0,0,0,-1,1,0))
#lrt15 <- glmLRT(fit, contrast = c(0,0,0,0,0,-1,0,1))
#lrt16 <- glmLRT(fit, contrast = c(0,0,0,0,0,0,-1,1))

table(decideTestsDGE(lrt01, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt02, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt03, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt04, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt05, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt06, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt07, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt08, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt09, p=0.05, adjust="BH"))
#table(decideTestsDGE(lrt10, p=0.05, adjust="BH"))
#table(decideTestsDGE(lrt11, p=0.05, adjust="BH"))
#table(decideTestsDGE(lrt12, p=0.05, adjust="BH"))
#table(decideTestsDGE(lrt13, p=0.05, adjust="BH"))
#table(decideTestsDGE(lrt14, p=0.05, adjust="BH"))
#table(decideTestsDGE(lrt15, p=0.05, adjust="BH"))
#table(decideTestsDGE(lrt16, p=0.05, adjust="BH"))

Toptable1 <- topTags(lrt01, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable2 <- topTags(lrt02, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable3 <- topTags(lrt03, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable4 <- topTags(lrt04, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable5 <- topTags(lrt05, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable6 <- topTags(lrt06, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable7 <- topTags(lrt07, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable8 <- topTags(lrt08, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
Toptable9 <- topTags(lrt09, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
#Toptable10 <- topTags(lrt10, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
#Toptable11 <- topTags(lrt11, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
#Toptable12 <- topTags(lrt12, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
#Toptable13 <- topTags(lrt13, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
#Toptable14 <- topTags(lrt14, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
#Toptable15 <- topTags(lrt15, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
#Toptable16 <- topTags(lrt16, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")