## Script 1: Variant Calling with MuTect2

#!/bin/bash

# Script to call variants with Mutect2

# This script will cd into a base folder and search for its subfolders

# cd into the base folder
cd /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2

# Find subfolders and store them into the $DIRS array variable
DIRS="$(find /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2/SRR1061196*/ -maxdepth 0 -type d)"

# Loop: for each folder in $DIRS
for j in ${DIRS[@]}
do
    cd $j

    # 1) Index the BAM files using SAMTOOLS
    samtools index *.bam

    # 2) Run MuTect2
    gatk Mutect2 -R /home/vgsouza/Vanessa_Doutorado/STAR/arquivos_de_entrada/hg38_genome.fa \
                 -L /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2/GetPileupSummaries/grep.af-only-gnomad.hg38.vcf \
                 -I bqsr_output.bam \
                 --germline-resource /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2/GetPileupSummaries/grep.af-only-gnomad.hg38.vcf \
                 --f1r2-tar-gz f1r2.tar.gz \
                 -O unfiltered.vcf \
                 -mbq 30

    # 3) Run LearnReadOrientationModel
    gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

    # 4) Run GetPileupSummaries to summarize read support for a set number of known variant sites
    gatk GetPileupSummaries -I bqsr_output.bam \
                             -V /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2/GetPileupSummaries/small_exac_common_3.hg38.vcf \
                             -L /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2/GetPileupSummaries/small_exac_common_3.hg38.vcf \
                             -O getpileupsummaries.table

    # 5) Estimate contamination with CalculateContamination
    gatk CalculateContamination -I getpileupsummaries.table \
                                 -tumor-segmentation segments.table \
                                 -O calculatecontamination.table

    # 6) FilterMutectCall
    gatk FilterMutectCalls -V unfiltered.vcf \
                            --tumor-segmentation segments.table \
                            --contamination-table calculatecontamination.table \
                            --ob-priors read-orientation-model.tar.gz \
                            -O filtered.vcf \
                            -R /home/vgsouza/Vanessa_Doutorado/STAR/arquivos_de_entrada/hg38_genome.fa

    # 7) Functional annotation with Funcotator
    gatk Funcotator -R /home/vgsouza/Vanessa_Doutorado/STAR/arquivos_de_entrada/hg38_genome.fa \
                    -V filtered.vcf \
                    -O outputFile \
                    --output-file-format VCF \
                    --data-sources-path /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2/funcotator_dataSources.v1.7.20200521s \
                    --ref-version hg38

    cd ../..

done



## Script 2: Separating SNPs and Indels with MuTect2


#!/bin/bash

# Script to separate SNPs and Indels using MuTect2

# This script will cd into a base folder and search for its subfolders

# cd into the base folder
cd /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2

# Find subfolders and store them into the $DIRS array variable
DIRS="$(find /home/vgsouza/Vanessa_Doutorado/Variant_Calling/MuTect2/SRR1061196*/ -maxdepth 0 -type d)"

# Loop: for each folder in $DIRS
for j in ${DIRS[@]}
do
    cd $j

    # Extract SNPs and store in mutect2.snp.vcf
    head -n 267 outputFile.pass.vcf > mutect2.snp.vcf && awk '{if ($8 ~ /SNP/) { print }}' outputFile.pass.vcf >> mutect2.snp.vcf

    # Extract Indels and store in mutect2.indel.vcf
    head -n 267 outputFile.pass.vcf > mutect2.indel.vcf && awk '{if ($8 ~ /INS/) { print }}' outputFile.pass.vcf >> mutect2.indel.vcf && awk '{if ($8 ~ /DEL/) { print }}' outputFile.pass.vcf >> mutect2.indel.vcf

    cd ../..

done
