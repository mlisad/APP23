####################################################################
# Author    : M. Dubbelaar
# Date      : 06-okt-2015
# File Name : loadingEnsemblData.R
# Purpose   : Creates a vector with all of the BioMart information 
#             of the specific gene.
####################################################################
getEnsembleData <- function(data) {
  data = M1
# Returns a list of the available Marts
listMarts()
# Defines the used Mart.
ensembl=useMart("ensembl")
# Returns the available datasets within ensembl. 
listDatasets(ensembl)
# Defines the dataset which will be used.
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)
#filters = listFilters(ensembl)

# These three lines describe the different data which can be used.
listDatasets(ensembl)[grep("musculus", listDatasets(ensembl)[,1]),]
attributes[grep("ensembl", attributes[,1]),]

# The getBM function retrieves the given attributes. 
BioM <- getBM(c("ensembl_gene_id", "gene_biotype","percentage_gc_content"), "", rownames(data), ensembl)
BioM2 <- getBM(c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"), "", rownames(data), ensembl)

# The gene names will be used to link the data to the data dataset.
# The entrez gene, the external gene name and the wiki gene description 
# can be put together with the matching ensembl gene id.
BioM <- BioM[match(rownames(data), BioM[,1]),]
BioM2 <- BioM2[match(rownames(data), BioM2[,1]),]
return (bioCom <- data.frame(BioM, BioM2[,2:4]))
}