####################################################################
# Author    : M. Dubbelaar
# Date      : 30-sept-2015
# File Name : linear_timepoints_deseq2.r
# Purpose   : The code within this file determines the genes who
#             play a role within the main effect of age, the main 
#             effect of genotype and the interaction effect of
#             genotype and age.
####################################################################
#                        Necessary functions                       #
####################################################################
calculateDds <- function(sampleTable, fileName1, fileName2, fileName3, dataType) {
  # This function calculates the dds.
  # DESeqDataSetFromHTSeqCount is used to store input values, intermediate 
  # calculations and results of the DE analysis.
  # The data will be filtered and normalized with the own made normalizeData function
  # The result from this function will be used within the DESeQ function
  # This function performs different statistical steps
  # The dds data will be used for the three other functions to determine the 
  # the main effect and the interaction effect.
  ddsHETSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~genotype*age)
  dds <- normalizeData(ddsHETSeq)
  dds <- DESeq(dds, betaPrior = F)
  
  calculateMainEffectGenotype(dds, fileName1)
  calculateMainEffectAge(dds, fileName2, dataType)
  calculateInteractionEffect(dds, fileName3, dataType)
}

calculateMainEffectGenotype <- function(dds, fileName){
  # The main effect of the genotype is for every age the same
  # This is the reason why this function is a small one.
  # The genotype of the WT is measured against the APP (HET) mice.
  # The results are transcribed into a txt file.
  mainEffectGenotype <- results(dds, name = "genotype_HET_vs_WT", alpha = 0.05)
  mainEffectGenotype <- getUniqueGenes(mainEffectGenotype)
  write.table(rownames(mainEffectGenotype), fileName, row.names = F, col.names=F, sep=", ", quote = F)
}

calculateMainEffectAge <- function(dds, fileName, dataType) {
  # The function calculateMainEffectAge needs to compare different ages.
  # The comparisons are made among all the ages, young mice (6-8 weeks and 12 months)
  # and old mice (12 months, 18 months and 24 months.)
  # The name of the comparison is used to get the genes and the unique
  # genes within these results are saved into a txt file.
  if (dataType != "old") {
    mainEffectAge1 <- results(dds, name = "age_6_vs_2", alpha = 0.05)
    mainEffectAge1 <- getUniqueGenes(mainEffectAge1)
    if (dataType == "all") {
      mainEffectAge2 <- results(dds, name =  "age_18_vs_2", alpha = 0.05)
      mainEffectAge2 <- getUniqueGenes(mainEffectAge2)
      mainEffectAge3 <- results(dds, name =  "age_24_vs_2",  alpha = 0.05)
      mainEffectAge3 <- getUniqueGenes(mainEffectAge3)
      mainEffectAge <- unique(c(row.names(mainEffectAge1), row.names(mainEffectAge2), row.names(mainEffectAge3)))
      write.table(mainEffectAge, fileName, row.names = F, col.names=F, sep=", ", quote = F)
    } else {
      write.table(rownames(mainEffectAge1), fileName, row.names = F, col.names=F, sep=", ", quote = F)
    }
  } else {
    mainEffectAge1 <- results(dds, name =  "age_18_vs_6", alpha = 0.05)
    mainEffectAge1 <- getUniqueGenes(mainEffectAge1)
    mainEffectAge2 <- results(dds, name =  "age_24_vs_6",  alpha = 0.05)
    mainEffectAge2 <- getUniqueGenes(mainEffectAge2)
    mainEffectAge <- unique(c(row.names(mainEffectAge1), row.names(mainEffectAge2)))
    write.table(mainEffectAge, fileName, row.names = F, col.names=F, sep=", ", quote = F)
  }
}

calculateInteractionEffect <- function(dds, fileName, dataType) {
  # The function calculateInteractionEffect needs to compare different ages.
  # The comparisons are made among all the ages, young mice (6-8 weeks and 12 months)
  # and old mice (12 months, 18 months and 24 months.)
  # This function is larger because of the different results in each age.
  # This data is also transcribed into a file.
  if (dataType != "young") {
    interaction1 <- results(dds, name="genotypeHET.age18", alpha=0.05)
    interaction1 <- getUniqueGenes(interaction1)
    interaction2 <- results(dds, name="genotypeHET.age24", alpha=0.05)
    interaction2 <- getUniqueGenes(interaction2)
    if (dataType == "all") {
      interaction3 <- results(dds, name="genotypeHET.age6", alpha=0.05)
      interaction3 <- getUniqueGenes(interaction1)
      interaction <- unique(c(row.names(interaction1), row.names(interaction2), row.names(interaction3)))
      write.table(interaction, fileName, row.names = F, col.names=F, sep=", ", quote = F)
    } else {
      interaction <- unique(c(row.names(interaction1), row.names(interaction2)))
      write.table(interaction, fileName, row.names = F, col.names=F, sep=", ", quote = F)
    }
  } else {
    interaction <- results(dds, name="genotypeHET.age6", alpha=0.05)
    interaction <- getUniqueGenes(interaction)
    write.table(rownames(interaction), fileName, row.names = F, col.names=F, sep=", ", quote = F)
  }
}

####################################################################
# The genotype is used to compared, it needs to be releved to make sure
# that the WT data is used first.
# The different ages of the mice are stored within the vector sampleAge
# After this the tables will be made with the values of the sampleGenotype
# and the sampleAge
sampleGenotype <- factor(targets$Genotype)
sampleGenotype <- relevel(sampleGenotype, "WT")
sampleAge <- factor(targets$Age)
####################################################################
#                             All ages                             #
####################################################################
# The data will be collected to create a table.
# This table is used by the function DESeqDataSetFromHTSeqCount to 
# create the right comparisons.
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, genotype=sampleGenotype, age=sampleAge)
calculateDds(sampleTable, "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/All_ages/main_genotype_result.txt", "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/All_ages/main_age_result.txt", 
                   "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/All_ages/interaction_result.txt", "all")
####################################################################
#                  12 months, 18 months and 24 months              #
####################################################################
sampleTable2 <- data.frame(sampleName=sampleFiles[7:24], fileName=sampleFiles[7:24], genotype=sampleGenotype[7:24], age=factor(targets$Age[7:24]))
calculateDds(sampleTable2, "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/6-18-24M_old_mice/main_genotype_result.txt", "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/6-18-24M_old_mice/main_age_result.txt", 
                            "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/6-18-24M_old_mice/interaction_result.txt", "old")
####################################################################
#                        6-8 weeks and 12 months                   #
####################################################################
sampleTable3 <- data.frame(sampleName=sampleFiles[1:12], fileName=sampleFiles[1:12], genotype=sampleGenotype[1:12], age=factor(sampleAge[1:12]))
calculateDds(sampleTable3, "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/2M-6M_old_mice/main_genotype_result.txt", "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/2M-6M_old_mice/main_age_result.txt", 
                           "/home/mdubbelaar/Desktop/APP23_results/DEseq2/Made_Documents/2M-6M_old_mice/interaction_result.txt", "young")
