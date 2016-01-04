####################################################################
# Author    : M. Dubbelaar
# Date      : 25-sept-2015
# File Name : plotColors.R
# Purpose   : Creating a color palette which can be used in multiple
#             files.
####################################################################

# A repetition will be made for each available target.
col_cell_age <- rep("black", dim(targets)[1])
# A contrast needs to be made for 2 different genotypes (APP and WT)
# and the age must differ in color as well, to determine the 
# differences between the genotype and the ages of the samples.
# The WT will be displayed with the use of a red color and the APP
# mice will be displayed with a blue color. The different 
# ages will differ from light (young mice) to dark (old mice).

# red color
col_cell_age[targets$Conditie ==  "WT.02"] = col=rainbow(265)[20]
col_cell_age[targets$Conditie ==  "WT.06"] = col=rainbow(265)[15]
col_cell_age[targets$Conditie ==  "WT.18"] = col=rainbow(265)[10]
col_cell_age[targets$Conditie ==  "WT.24"] = col=rainbow(265)[1]
# Blue color
col_cell_age[targets$Conditie ==  "HET.02"] = col=rainbow(265)[145]
col_cell_age[targets$Conditie ==  "HET.06"] = col=rainbow(265)[150]
col_cell_age[targets$Conditie ==  "HET.18"] = col=rainbow(265)[155]
col_cell_age[targets$Conditie ==  "HET.24"] = col=rainbow(265)[160]
