#!/bin/bash

## Script to annotate using Annovar and VarScan2

# Navigate to the base folder
cd /home/vgsouza/Vanessa_Doutorado/Variant_Calling/VarScan_results

# Find subfolders and store them in the $DIRS array variable
DIRS=$(find /home/vgsouza/Vanessa_Doutorado/Variant_Calling/VarScan_results/SRR1061196*/ -maxdepth 0 -type d)

# Loop through each folder in $DIRS
for j in ${DIRS[@]}
do
    cd $j

    ## Annotate SNP files with Annovar
    /home/vgsouza/Vanessa_Doutorado/annovar/annovar/./table_annovar.pl varScan.snp.annovar /home/vgsouza/Vanessa_Doutorado/annovar/annovar/humandb/ -build hg38 -out ex4 --outfile annovar_snp -remove -protocol knownGene,dbnsfp30a,cosmic92,clinvar_20210501,avsnp150 -operation g,f,f,f,f -nastring . -polish

    ## Annotate INDEL files with Annovar
    /home/vgsouza/Vanessa_Doutorado/annovar/annovar/./table_annovar.pl varScan.indel.annovar /home/vgsouza/Vanessa_Doutorado/annovar/annovar/humandb/ -build hg38 -out ex4 --outfile annovar_indel -remove -protocol knownGene,dbnsfp30a,cosmic92,clinvar_20210501,avsnp150 -operation g,f,f,f,f -nastring . -polish

    cd ..
done
