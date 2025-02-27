#!/bin/bash

## Script for annotating with Annovar Mutect2 results

# This script will change directory to a base folder, search for its subfolders

# Change directory to base folder
cd /home/vgsouza/Vanessa_Doutorado/Variant_Calling/Mutect2_results

# Find subfolders and store them into $DIRS array variable
DIRS=$(find /home/vgsouza/Vanessa_Doutorado/Variant_Calling/Mutect2_results/SRR1061196*/ -maxdepth 0 -type d)

# Loop: for each folder in $DIRS
for j in ${DIRS[@]}
do
    cd $j

    ## Annovar for SNP files
    /home/vgsouza/Vanessa_Doutorado/annovar/annovar/./table_annovar.pl mutect.snp.annovar /home/vgsouza/Vanessa_Doutorado/annovar/annovar/humandb/ -build hg38 -out ex4 --outfile annovar -remove -protocol knownGene,dbnsfp30a,cosmic92,clinvar_20210501,avsnp150 -operation g,f,f,f,f -nastring . -
