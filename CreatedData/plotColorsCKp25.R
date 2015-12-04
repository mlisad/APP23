####################################################################
# Author    : M. Dubbelaar
# Date      : 4-dec-2015
# File Name : plotColorsCKp25.R
# Purpose   : Creating a color palette which can be used in multiple
#             files.
####################################################################

# A repetition will be made for each available target.
col_cell_age <- rep("black", dim(targets)[1])
col_cell_age[targets$Conditie ==  "WT_02w"] = col=rainbow(265)[20]
col_cell_age[targets$Conditie ==  "WT_06w"] = col=rainbow(265)[1]
# Blue color
col_cell_age[targets$Conditie ==  "CKp25_02W"] = col=rainbow(265)[145]
col_cell_age[targets$Conditie ==  "CKp25_06W"] = col=rainbow(265)[160]