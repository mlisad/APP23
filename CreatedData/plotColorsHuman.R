####################################################################
# Author    : M. Dubbelaar
# Date      : 04-dec-2015
# File Name : plotColorsHuman.R
# Purpose   : Creating a color palette which can be used in multiple
#             files.
####################################################################

# A repetition will be made for each available target.
col_cell_age <- rep("black", dim(targets)[1])
col_cell_age[targets$Genotype ==  "LOAD"] = col=rainbow(265)[20]
col_cell_age[targets$Genotype ==  "CTRL"] = col=rainbow(265)[145]