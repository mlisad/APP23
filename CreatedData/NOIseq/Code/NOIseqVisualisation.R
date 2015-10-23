####################################################################
# Author    : M. Dubbelaar
# Date      : 05-okt-2015
# File Name : NOIseqVisualisation.R
# Purpose   : To visualise the APP data with the package NOIseq.
# Used Files: loadingAppFile.R
#             loadingEnsemblDataNew.R
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/Desktop/Onderzoek-APP23_RNASEQ/CreatedData/NOIseq/")
#biocLite("NOISeq")
#biocLite("biomaRt")
library("NOISeq")
library("biomaRt")
####################################################################
#               Reading of the data and the targets                #
####################################################################
source("../loadingAppFile.R")
# The data found in loadingAppFile needs some adjustment for further
# progress. The data needs to be filtered of the genenames within 
# the column and needs to be set as the column name.
M1 <- APP23_data[,2:25]
row.names(M1) <- APP23_data[,1]
####################################################################
#                      Loading Ensembl data                        #
####################################################################
source("Code/loadingEnsemblDataNEW.R")

createFunction <- function(data, pathPlot, pathDoc, type){
  BioM <- getEnsembleData(data)
####################################################################
#                    Creating filtered dataset                     #
####################################################################
# The following data needs to be defined because of the 
# rownames for each vector. When the vector is not created
# and the data which can be found in the data frame is defined
# directly, the first step will fail.
myLength <- data.frame(BioM[,6]-BioM[,5], row.names = rownames(data))
myGc <- data.frame(BioM[,3], row.names = rownames(data))
myBiotypes <- data.frame(BioM[,2], row.names = rownames(data))
myChromosomes <- data.frame(BioM[,4:6], row.names = rownames(data))

if (type == "All") {
  Genotype <- targets$Genotype
  myfactors <- data.frame(Genotype, targets$Conditie, targets$Samples, targets$Age)
} else if (type == "Old"){
  Genotype <- targets$Genotype[7:24]
  myfactors <- data.frame(Genotype, targets$Conditie[7:24], targets$Samples[7:24], targets$Age[7:24])
} else if (type == "Young") {
  Genotype <- targets$Genotype[1:12]
  myfactors <- data.frame(Genotype, targets$Conditie[1:12], targets$Samples[1:12], targets$Age[1:12])
}
myFilt <- filtered.data(data, factor = myfactors$targets.Conditie, method = 3, norm = F, cpm=2, p.adj = "BH")
myData <- readData(data=myFilt, length=myLength, gc=myGc, biotype = myBiotypes,
                              chromosome = myChromosomes, factors = myfactors)
####################################################################
#                       Quality control                            #
####################################################################
# Within these lines of code the data can be visualized on different
# features. Some exaples are the count distribution per biotype.
# Where the biotype of each sample will be displayed with the number of
# the detected features displayed in the top of the plot.
myCountsBio <- dat(myData, factor = NULL, type = "countsbio")
pdf(paste(pathPlot, "biotypeCounts.pdf", sep=''))
if (type == "All") {
  explo.plot(myCountsBio, toplot = 1, samples = c(1:24), plottype = "boxplot")
  # The plot below shows the different biotypes within one sample.
  explo.plot(myCountsBio, toplot = 1, samples = 19, plottype = "boxplot")
} else if (type == "Old"){
  explo.plot(myCountsBio, toplot = 1, samples = c(1:18), plottype = "boxplot")
  # The plot below shows the different biotypes within one sample.
  explo.plot(myCountsBio, toplot = 1, samples = 18, plottype = "boxplot")
} else if (type == "Young") {
  explo.plot(myCountsBio, toplot = 1, samples = c(1:12), plottype = "boxplot")
  # The plot below shows the different biotypes within one sample.
  explo.plot(myCountsBio, toplot = 1, samples = 12, plottype = "boxplot")
}
# The plot below shows the amount of the group protein_coding 
# for each sample.
explo.plot(myCountsBio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
# The plot below shows how high the count per milion in each sample is.
# this cpm are divided in different groups. 
explo.plot(myCountsBio, toplot = 1, samples = NULL, plottype = "barplot")
dev.off()
# The saturation shows the amound of features that are detected in the genome 
# that exceed the number k. 
mySaturation <- dat(myData, k = 10, ndepth = 7, type = "saturation")
pdf(paste(pathPlot, "saturation.pdf", sep=''))
explo.plot(mySaturation, toplot = 1, samples = c(7,10), yleftlim = NULL, yrightlim = NULL)
dev.off()
# The lengthbias is used to describe to relation between the length and 
# the expression values. When the p value is significant and the
# R2 value is above the 70% The expression will depend on the length.
# The curve shows the type of dependence.
myLengthBias <- dat(myData, factor = "Genotype", type = "lengthbias")
pdf(paste(pathPlot, "lengthBias.pdf", sep=''))
explo.plot(myLengthBias, samples = NULL, toplot = "global")
dev.off()
# The GCbias describes the relationship betwen the GC content and the 
# expression values. Again when the p value is significant and the R2
# value is higher than 70% the expression will depend on the the GC content
myGCbias <- dat(myData, factor = "Genotype", type = "GCbias")
pdf(paste(pathPlot, "GCBias.pdf", sep=''))
explo.plot(myGCbias, samples = NULL, toplot = "global")
dev.off()
# The CD plot shows a RNA composition plot.
myCD <- dat(myData, type = "cd", norm = F, refColumn = 1)
pdf(paste(pathPlot, "cd.pdf", sep=''))
if (type == "All") {
  explo.plot(myCD, samples = 1:12)     
  explo.plot(myCD, samples = 13:24)
} else if (type == "Old") {
  explo.plot(myCD, samples = 1:9)
  explo.plot(myCD, samples = 10:18)
} else if (type == "Young") {
  explo.plot(myCD, samples = 1:12)
}
dev.off()
####################################################################
#                            Normalization                         #
####################################################################
# The data can be normalized with the following techniques
# upperquartile (uqua), RPKM and Trimmed Mean of M values (TMM) 
myRPKM <- rpkm(assayData(myData)$exprs, lc = 1, k = 0)
myUQUA <- uqua(assayData(myData)$exprs, lc = 0.5, k = 0)
myTMM <- tmm(assayData(myData)$exprs, long = 1000, lc = 0)
####################################################################
#                               Results                            #
####################################################################
# The estimate of the probability distribution for M and D are made.
# M and D values for every pair of replicates with the same condition
# are computed. noiSeqBio is optimized to use on biological replicates.
myNoiSeqBio <- noiseqbio(myData, k = 0.5, norm = "rpkm", factor = "Genotype",
                         lc = 0, r = 15, plot = F, filter = 2)
# The differential expressed features will be obtained
# with the use of the function degenes.
# The genes will be grouped when they are up or down regulated.
myNoiSeq.degAll <- degenes(myNoiSeqBio, q =  0.7, M = NULL)
myNoiSeq.degUp <- degenes(myNoiSeqBio, q =  0.7, M = "up")
myNoiSeq.degDown <- degenes(myNoiSeqBio, q =  0.7, M = "down")
pdf(paste(pathPlot, "DEPlots.pdf", sep=''))
# The expression plot is shown
DE.plot(myNoiSeqBio, q = 0.7, graphic = "expr", log.scale = T)
# The MD plot is shown, the log fold change (M) and the absolute
# value of the difference in expression between the conditions (D)
# are plotted against each other.
DE.plot(myNoiSeqBio, q = 0.7, graphic = "MD")
# The Manhattan plots are made for the chromosomes.
# This plot shows the up of down regulation of genes.
DE.plot(myNoiSeqBio, chromosomes = c(1,4, 10), log.scale = T, join = F, q = 0.5, graphic = "chrom")
DE.plot(myNoiSeqBio, chromosomes = c(11, 15, 16), log.scale = T, join = F, q = 0.5, graphic = "chrom")
DE.plot(myNoiSeqBio, chromosomes = c(18, 19), log.scale = T, join = F, q = 0.5, graphic = "chrom")
# Regulation on both genotypes on the chromosome.
DE.plot(myNoiSeqBio, chromosomes = c(2,5), log.scale = T, join = F, q = 0.5, graphic = "chrom")
DE.plot(myNoiSeqBio, chromosomes = c(7,9, 17), log.scale = T, join = F, q = 0.5, graphic = "chrom")

# The distribution of differentially expressed features of each chromosome
# and biotypes are shown.
DE.plot(myNoiSeqBio, chromosomes = c(1:19), q = 0.7, graphic = "distr")
dev.off()

write.table(rownames(myNoiSeq.degAll), paste(pathDoc, "DEPlotsAllReg.txt", sep=""), row.names = F, col.names=F, eol=",\n", quote = F)
write.table(rownames(myNoiSeq.degUp), paste(pathDoc, "DEPlotsUpReg.txt",  sep=""), row.names = F, col.names=F, eol=",\n", quote = F)
write.table(rownames(myNoiSeq.degDown), paste(pathDoc, "DEPlotsDownReg.txt",  sep=""), row.names = F, col.names=F, eol=",\n", quote = F)
}
####################################################################
createFunction(M1, "Plots/All_mice/", "Made_Documents/All_ages/", "All")
createFunction(M1[,7:24], "Plots/Old_mice/", "Made_Documents/12-18-24M_old_mice/", "Old")
createFunction(M1[,1:12], "Plots/Young_mice/", "Made_Documents/6_8M-12M_old_mice/","Young")