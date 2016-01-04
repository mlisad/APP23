####################################################################
# Author    : M. Dubbelaar
# Date      : 30-sept-2015
# File Name : linear_timepoints_deseq2.r
# Purpose   : The code within this file determines the genes who
#             play a role within the main effect of age, the main 
#             effect of genotype and the interaction effect of
#             genotype and age.
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
calculateDds(sampleTable, paste(resultPathway, "Made_Documents/All_ages/main_genotype_result.txt", sep=""), paste(resultPathway, "Made_Documents/All_ages/main_age_result.txt", sep=""), 
                   paste(resultPathway, "Made_Documents/All_ages/interaction_result.txt", sep=""), "all")
####################################################################
#                  12 months, 18 months and 24 months              #
####################################################################
sampleTable2 <- data.frame(sampleName=sampleFiles[7:24], fileName=sampleFiles[7:24], genotype=sampleGenotype[7:24], age=factor(targets$Age[7:24]))
calculateDds(sampleTable2, paste(resultPathway, "Made_Documents/6-18-24M_old_mice/main_genotype_result.txt", sep=""), paste(resultPathway, "Made_Documents/6-18-24M_old_mice/main_age_result.txt", sep=""), 
                            paste(resultPathway, "Made_Documents/6-18-24M_old_mice/interaction_result.txt", sep=""), "old")
####################################################################
#                        6-8 weeks and 12 months                   #
####################################################################
sampleTable3 <- data.frame(sampleName=sampleFiles[1:12], fileName=sampleFiles[1:12], genotype=sampleGenotype[1:12], age=factor(sampleAge[1:12]))
calculateDds(sampleTable3, paste(resultPathway, "Made_Documents/2M-6M_old_mice/main_genotype_result.txt", sep=""), paste(resultPathway, "Made_Documents/2M-6M_old_mice/main_age_result.txt", sep=""), 
                           paste(resultPathway, "Made_Documents/2M-6M_old_mice/interaction_result.txt", sep=""), "young")
