#!/bin/bash
#
# The first step is to download SRAToolkit
# http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
#
# This scipt is used to create a fastq file that is ready for the 
# fastQC. The first step is to load the SRA file, this is the input that
# must be given while calling the file. This data must be provided as a 
# string (example: "SRR1127223 SRR1127224 SRR1127225")
# These .sra files will be transformed into fastq files later on.

############################################################
#					      Pathways						   #
############################################################
# The pathways can be changed to the pathway of your own.
sraPackagePathway="/home/mdubbelaar/Desktop/AligningAndQC/sratoolkit.2.5.6-ubuntu64/bin"
sraFilesPathway="/home/mdubbelaar/ncbi/public/sra"
fastQFilesPathway="/media/mdubbelaar/BD_4T/FastQ/"
############################################################
#							Code						   #
############################################################
# For each element in the given string a prefetch will be done
# This function makes sure that the .sra files are collected.  
for sraNr in $1;
do
	echo "Getting the SRA files; $sraNr"
	$sraPackagePathway/prefetch $sraNr
# Takes the 6th element of a spliced fastq string, which is the filename.
done

# For each element in the sra folder.
# The item will be transcribed into a fastQ file in the given directory.
for sraFile in $sraFilesPathway/*.sra;
do 	
	echo "Writing $sraFile as a fastQ file"
	$sraPackagePathway/fastq-dump $sraFile --outdir $fastQFilesPathway
# Takes the 6th element of a spliced fastq string, which is the filename.
done	


