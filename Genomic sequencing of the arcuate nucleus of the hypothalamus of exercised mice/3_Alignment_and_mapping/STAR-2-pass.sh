#!/bin/bash
## script mapping using STAR-2-pass

## 1) Indexing genome with annotations
## Indexing for maximum read lengh 100 bp

/home/vgsouza/STAR/STAR-2.7.11b/source/./STAR --runMode genomeGenerate --genomeDir /home/vgsouza/VanessaPD/reference/ --genomeFastaFiles /home/vgsouza/VanessaPD/reference/GCF_000001635.27_GRCm39_genomic.fna --sjdbGTFfile /home/vgsouza/VanessaPD/reference/GCF_000001635.27_GRCm39_genomic.gtf --runThreadN 35 --sjdbOverhang 99


## 2) 1-pass mapping with indexed genome

# Go to the "samples" directory
cd /home/vgsouza/VanessaPD/samples

# Loop through directories 1 to 20
for i in {1..20}; do
    # Checks if the directory exists
    if [ -d "$i" ]; then
        # Enter the "Thinned" directory
        cd "$i/Thinned" || exit
        
        # run star 1-pass mapping with indexed genome
        /home/vgsouza/STAR/STAR-2.7.11b/source/./STAR --genomeDir /home/vgsouza/VanessaPD/reference/ --readFilesIn seqyclean_PE1.fastq seqyclean_PE2.fastq --runThreadN 40 --outSAMunmapped Within --outFileNamePrefix ./Prefix
        
        # Back to the "samples" directory
        cd ../../
    fi
done


### NOTE ###
## The same command has been run for multiple samples in the for loop, therefore, it will generate SJ.out.tab file for each sample.
## Next, I have copied SJ.out.tab files of all the samples into the a single folder "SJ_out"

mkdir /home/vgsouza/VanessaPD/SJ_out

cd /home/vgsouza/VanessaPD/samples

# Loop through directories 1 to 20
for i in {1..20}; do
    # Checks if the directory exists
    if [ -d "$i" ]; then
        # Enter the "Thinned" directory
        cd "$i/Thinned" || exit
        
        # cp PrefixSJ.out.tab
	cp *PrefixSJ.out.tab /home/vgsouza/VanessaPD/SJ_out
        
        # Back to the "samples" directory
        cd ../../
    fi
done



## 3) Indexing genome with annotations and SJ.out.tab files
## Note: Again indexing for maximum read lengh 90 bp.

/home/vgsouza/STAR/STAR-2.7.11b/source/./STAR --runMode genomeGenerate --genomeDir /home/vgsouza/VanessaPD/reference/SJ_Index/ --genomeFastaFiles /home/vgsouza/VanessaPD/reference/SJ_Index/GCF_000001635.27_GRCm39_genomic.fna --sjdbGTFfile /home/vgsouza/VanessaPD/reference/SJ_Index/GCF_000001635.27_GRCm39_genomic.gtf --runThreadN 40 --sjdbOverhang 99 --sjdbFileChrStartEnd /home/vgsouza/VanessaPD/SJ_out/*PrefixSJ.out.tab


## 4) 2-pass mapping with new indexed genome with annotations and SJ.out.tab files

#!/bin/bash
cd /home/vgsouza/VanessaPD/samples

# Loop through directories 1 to 20
for i in {1..20}; do
    # Checks if the directory exists
    if [ -d "$i" ]; then
        # Enter the "Thinned" directory
        cd "$i/Thinned" || exit
        
        # run star 2-pass mapping with indexed genome
        /home/vgsouza/STAR/STAR-2.7.11b/source/./STAR --genomeDir /home/vgsouza/VanessaPD/reference/SJ_Index/ --readFilesIn seqyclean_PE1.fastq seqyclean_PE2.fastq --runThreadN 40 --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix ./Prefix2
        
        # Back to the "samples" directory
        cd ../../
    fi
done