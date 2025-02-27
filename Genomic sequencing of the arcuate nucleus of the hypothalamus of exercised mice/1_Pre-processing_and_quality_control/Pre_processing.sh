#!/bin/bash

# This script will cd into a base folder, search for its subfolders

# cd into base folder
cd /home/vgsouza/VanessaPD/samples

# find subfolders and store them into $DIRS array variable
DIRS="$(find ./*/ -maxdepth 0 -type d)"

# loop: for each folder in $DIRS
for j in ${DIRS[@]}
do
    (
    cd $j

    #run fastqc
    nohup fastqc *.fastq.gz &

    #run thinning with seqyclean
    nohup seqyclean -v /home/vgsouza/Vanessa_Doutorado/UniVec -qual 24 24 -minlen 65 -1 *R1_001.fastq.gz -2 *R2_001.fastq.gz -o Thinned/seqyclean -t 20 &

    # run fastqc
    nohup fastqc Thinned/*.fastq &

    ## move the fastqc.zip BEFORE to cleanning to a folder called before_cleaning
    nohup mv *fastqc.zip /home/vgsouza/VanessaPD/quality_data/before_cleaning/*_fastqc.zip &

    ## run multiqc
    nohup multiqc /home/vgsouza/VanessaPD/quality_data/before_cleaning/ &

    ## move the fastqc.zip AFTER to cleanning to a folder called after_cleaning
    cd /home/vgsouza/VanessaPD/samples
    for i in {1..20}; do
    # Checks if the directory exists
    if [ -d "$i" ]; then
        # Enter the "Thinned" directory
        cd "$i/Thinned" || exit
        
        # Move *fastqc.zip files to specific directory
        mv *fastqc.zip /home/vgsouza/VanessaPD/quality_data/after_cleaning
        
        # Returns to the "samples" directory
        cd ../../
    fi
done

    ## run multiqc
    nohup multiqc /home/vgsouza/VanessaPD/quality_data/after_cleaning/ &

    cd ..
    ) &
done

wait


