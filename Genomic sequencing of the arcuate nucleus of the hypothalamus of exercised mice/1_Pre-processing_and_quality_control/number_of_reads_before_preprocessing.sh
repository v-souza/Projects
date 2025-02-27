#!/bin/bash

# Define o diretório raiz
root_dir="/home/vgsouza/VanessaPD/samples"

# Define o arquivo onde o resultado será salvo
output_file="resultados.txt"

# Remove o arquivo de resultados se existir
rm -f "$output_file"

# Loop através dos diretórios numerados de 1 a 20
for dir in {1..20}; do
    # Verifica se o diretório existe
    if [ -d "$root_dir/$dir" ]; then
        # Loop através dos arquivos *fastq.gz
        for file in "$root_dir/$dir"/*fastq.gz; do
            # Calcula o número de linhas dividido por 4
            count=$(zcat "$file" | wc -l)
            result=$((count / 4))
            
            # Obtém o nome do arquivo sem o caminho completo
            filename=$(basename "$file")
            
            # Salva o nome do arquivo e o resultado no arquivo de saída
            echo -e "$filename\t$result" >> "$output_file"
        done
    fi
done
