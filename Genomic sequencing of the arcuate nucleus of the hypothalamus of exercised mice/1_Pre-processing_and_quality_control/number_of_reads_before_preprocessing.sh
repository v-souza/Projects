#!/bin/bash

# Define o diret�rio raiz
root_dir="/home/vgsouza/VanessaPD/samples"

# Define o arquivo onde o resultado ser� salvo
output_file="resultados.txt"

# Remove o arquivo de resultados se existir
rm -f "$output_file"

# Loop atrav�s dos diret�rios numerados de 1 a 20
for dir in {1..20}; do
    # Verifica se o diret�rio existe
    if [ -d "$root_dir/$dir" ]; then
        # Loop atrav�s dos arquivos *fastq.gz
        for file in "$root_dir/$dir"/*fastq.gz; do
            # Calcula o n�mero de linhas dividido por 4
            count=$(zcat "$file" | wc -l)
            result=$((count / 4))
            
            # Obt�m o nome do arquivo sem o caminho completo
            filename=$(basename "$file")
            
            # Salva o nome do arquivo e o resultado no arquivo de sa�da
            echo -e "$filename\t$result" >> "$output_file"
        done
    fi
done
