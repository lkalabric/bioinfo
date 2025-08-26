#!/bin/bash

# Uso do bepipred
# Link: https://github.com/UberClifford/BepiPred3.0-Predictor/tree/main

# BepiPred3.0 predicts B-cell epitopes from ESM-2 encodings of proteins sequences.

# Instalação
conda create -n bepipred3 python=3.8.8
conda activate bepipred3
conda install pip
pip3 install -r requirements.txt

# Define o nome do arquivo de entrada e de saÃ­da
#input_file="hdv-1_bepipred.txt"
input_file=$1
#output_file="hdv-1_bepipred.sam"
output_file="${input_file%.*}.sam"

# Verifica se o arquivo de entrada existe
if [[ ! -f "$input_file" ]]; then
    echo "Erro: O arquivo '$input_file' não foi encontrado."
    exit 1
fi

# Usa awk para processar o arquivo e gerar a saÃ­da no formato .sam.
# Redireciona a saÃ­da completa do awk para o arquivo de saÃ­da.
# `last_col` Ã© a Ãºltima coluna, o 9Âº campo.
# `start_col` Ã© a 4Âª coluna.
# `end_col` Ã© a 5Âª coluna.
awk '
BEGIN {
    # Define o delimitador de campo como espaÃ§o
    FS = " "  

    # Imprime o cabeÃ§alho do arquivo SAM
    print "@HD\tVN:1.6\tSO:unknown"
    print "@SQ\tSN:AEV40786.1\tLN:106" # Exemplo, o comprimento deve ser ajustado
    print "@PG\tID:bepipred-parser\tPN:bepipred-parser\tVN:1.0"
    
    in_block = 0  # Flag para indicar se estamos dentro de um bloco de "E"
    start_pos = 0 # VariÃ¡vel para armazenar a posiÃ§Ã£o de inÃ­cio
    count = 0     # Contador para o nÃºmero de "E"s consecutivos
}

# Pula as linhas de cabeÃ§alho
/^#/ {next}
/^-/ {next}

{
    last_col = $9
    start_pos_current = $4
    end_pos_current = $5

    if (last_col == "E") {
        # Se a linha atual tiver "E"
        if (in_block == 0) {
            # Inicia um novo bloco e reinicia a contagem
            in_block = 1
            start_pos = start_pos_current
            count = 1
        } else {
            # Incrementa a contagem se jÃ¡ estiver em um bloco
            count++
        }
        # Atualiza a posiÃ§Ã£o de fim para a linha atual
        end_pos = end_pos_current
    } else {
        # Se a linha atual nÃ£o tiver "E"
        if (in_block == 1) {
            # O bloco terminou.
            # Verifica se a contagem atende ao requisito mÃ­nimo de 8
            if (count >= 8) {
                # Se sim, imprime a linha no formato SAM
                # Colunas SAM: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
                print "AEV40786.1\t0\thdv-1\t" start_pos "\t255\t*\t*\t0\t0\t*\t*"
            }
            # Reseta as variÃ¡veis
            in_block = 0
            start_pos = 0
            count = 0
        }
    }
}

END {
    # Se o arquivo terminar dentro de um bloco, verifica e imprime o Ãºltimo bloco
    if (in_block == 1 && count >= 8) {
        print "AEV40786.1\t0\thdv-1\t" start_pos "\t255\t*\t*\t0\t0\t*\t*"
    }
}
' "$input_file" > "$output_file"

echo "Arquivo '$output_file' gerado com sucesso!"
