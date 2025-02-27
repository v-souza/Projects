#!/bin/bash

# This script performs variant calling, filtering, and conversion to ANNOVAR format

# Define base folder
BASE_FOLDER="/home/vgsouza/Vanessa_Doutorado/Variant_Calling/VarScan2"

# Function to process each sample folder
process_sample() {
    local sample_folder="$1"
    
    # Move into sample folder
    cd "$sample_folder"

    # Sort BAM file
    samtools sort bqsr_output.bam > bqsr_output.sorted.bam
    
    # Index BAM file
    samtools index bqsr_output.sorted.bam
    
    # Generate pileup
    samtools mpileup --min-MQ 1 --no-BAQ --fasta-ref /home/vgsouza/Vanessa_Doutorado/STAR/arquivos_de_entrada/hg38_genome.fa bqsr_output.sorted.bam > myData.pileup
    
    # Call SNPs
    varscan mpileup2snp myData.pileup --min-coverage 10 --min-var-freq 0.20 --p-value 0.05 --output-vcf 1 > sample.varScan.snp
    
    # Call Indels
    varscan mpileup2indel myData.pileup --min-coverage 10 --min-var-freq 0.10 --p-value 0.10 --output-vcf 1 > sample.varScan.indel
    
    # Filter SNPs near Indels
    varscan filter sample.varScan.snp --indel-file sample.varScan.indel --output-file sample.varScan.snp.filter
    
    # Filter Indels for higher confidence
    varscan filter sample.varScan.indel --min-reads2 4 --min-var-freq 0.15 --p-value 0.05 --output-file sample.varScan.indel.filter

    # Filter SNPs and Indels passing false positive filter
    head -n 36 sample.varScan.snp.fpfilter > sample.varScan.snp.fpfilter.pass.vcf && awk '{ if($7 == "PASS") { print }}' sample.varScan.snp.fpfilter >> sample.varScan.snp.fpfilter.pass.vcf
    head -n 36 sample.varScan.indel.fpfilter > sample.varScan.indel.fpfilter.pass.vcf && awk '{ if($7 == "PASS") { print }}' sample.varScan.indel.fpfilter >> sample.varScan.indel.fpfilter.pass.vcf
    
    # Convert to ANNOVAR format
    /home/vgsouza/Vanessa_Doutorado/annovar/annovar/convert2annovar.pl -format vcf4 -includeinfo sample.varScan.snp.fpfilter.pass.vcf > sample.varScan.snp.fpfilter.pass.vcf.annovar
    /home/vgsouza/Vanessa_Doutorado/annovar/annovar/convert2annovar.pl -format vcf4 -includeinfo sample.varScan.indel.fpfilter.pass.vcf > sample.varScan.indel.fpfilter.pass.vcf.annovar
    
    # Run fpfilter
    perl /home/vgsouza/Vanessa_Doutorado/fpfilter-tool-master/fpfilter.pl \
      --vcf-file sample.varScan.snp.filter.vcf \
      --bam-file bqsr_output.sorted.bam \
      --bam-index bqsr_output.sorted.bam.bai \
      --reference /home/vgsouza/Vanessa_Doutorado/STAR/arquivos_de_entrada/hg38_genome.fa \
      --sample Sample1 \
      --output sample.varScan.snp.fpfilter

    perl /home/vgsouza/Vanessa_Doutorado/fpfilter-tool-master/fpfilter.pl \
      --vcf-file sample.varScan.indel.filter.vcf \
      --bam-file bqsr_output.sorted.bam \
      --bam-index bqsr_output.sorted.bam.bai \
      --reference /home/vgsouza/Vanessa_Doutorado/STAR/arquivos_de_entrada/hg38_genome.fa \
      --sample Sample1 \
      --output sample.varScan.indel.fpfilter

    # Move back to parent folder
    cd "$BASE_FOLDER"
}

# Main loop: process each sample folder
for sample_dir in $(find "$BASE_FOLDER" -type d -name "SRR1061196*"); do
    process_sample "$sample_dir"
done
