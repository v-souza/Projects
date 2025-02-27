#!/bin/bash

# Sets the root directory
root_dir="/home/vgsouza/VanessaPD/samples"

# Defines the file where the result will be saved
output_file="after_preprocessing.txt"

# Remove the results file if it exists
rm -f "$output_file"

# Loop through directories numbered 1 to 20
for dir in {1..20}; do
    # Checks if the directory exists
    if [ -d "$root_dir/$dir/Thinned" ]; then
        # Loop through *fastq files in Thinned
        for file in "$root_dir/$dir/Thinned"/*fastq; do
            # calculates the number of lines divided by 4
            count=$(cat "$file" | wc -l)
            result=$((count / 4))
            
            # Gets the file name without the full path
            filename=$(basename "$file")
            
            # Saves the file name and result to the output file
            echo -e "$filename\t$result" >> "$output_file"
        done
    fi
done
