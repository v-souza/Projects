#!/bin/bash
## script Counting reads using htseq-count


# Go to the "samples" directory
cd /home/vgsouza/VanessaPD/samples

# Loop through directories 1 to 20
for i in {1..20}; do
    # Checks if the directory exists
    if [ -d "$i" ]; then
        # Enter the "Thinned" directory
        cd "$i/Thinned" || exit
     
## create index with samtools
samtools index *.bam

## run htseq-count
htseq-count -m intersection-nonempty -i gene_id -r pos -s no -f bam Prefix2Aligned.sortedByCoord.out.bam /home/vgsouza/VanessaPD/reference/SJ_Index/GCF_000001635.27_GRCm39_genomic.gtf > htseq-count.txt


# Back to the "samples" directory
        cd ../../
    fi
done