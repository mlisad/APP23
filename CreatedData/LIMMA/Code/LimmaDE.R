####################################################################
# Author    : M. Dubbelaar
# Date      : 22-sept-2015
# File Name : LimmaDE.R
# Purpose   : Calculates the values of the estimates which where created by limma.
####################################################################
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