#!/bin/bash

# Este script converte um arquivo CSV de mapeamento simplificado para o formato SAM.
# Ele assume que o seu CSV tem um cabeçalho e as seguintes colunas (nesta ordem):
# read_id, ref_id, posicao, sequencia
# Substitua 'seu_arquivo.csv' pelo nome do seu arquivo de entrada.

# Nome do arquivo de entrada CSV
INPUT_CSV=$1

# Nome do arquivo de saída SAM
OUTPUT_SAM="${INPUT_CSV%.*}.sam"

# --- 1. Criar o cabeçalho do arquivo SAM ---
# O cabeçalho é obrigatório para um arquivo SAM válido.
# Aqui, criamos um cabeçalho simples.
SAM_HEADER="@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref1\tLN:1000\n@PG\tID:my_script\tPN:my_script\tVN:1.0"

echo -e "$SAM_HEADER" > "$OUTPUT_SAM"

# --- 2. Processar o arquivo CSV linha por linha ---
# Usamos 'tail -n +2' para pular o cabeçalho do CSV.
# O comando 'while read' lê cada linha do CSV.
# O 'IFS=,' define o separador de campo como vírgula.
tail -n +2 "$INPUT_CSV" | while IFS=',' read -r read_id ref_id posicao sequencia
do
    # --- 3. Preencher os campos do registro SAM ---
    # Campos fixos ou calculados.
    flag=0
    mapq=20
    ref_id="ref1"
    # A string CIGAR é o comprimento da sequência seguido de 'M' (Match).
    cigar="${#sequencia}M"
    rnext="*"
    pnext=0
    tlen=0
    # A qualidade (QUAL) é definida como '*' pois não está disponível no CSV.
    qual="*"

    posicao_corrigida=$((posicao-1))
    # --- 4. Construir o registro SAM (separado por TABs) ---
    # sam_record="${read_id}\t${flag}\t${ref_id}\t${posicao}\t${mapq}\t${cigar}\t${rnext}\t${pnext}\t${tlen}\t${sequencia}\t${qual}"
    sam_record="${read_id}\t${flag}\t${ref_id}\t${posicao_corrigida}\t${mapq}\t${cigar}\t${rnext}\t${pnext}\t${tlen}\t${sequencia}\t${qual}"
    # --- 5. Escrever o registro no arquivo de saída ---
    echo -e "$sam_record" >> "$OUTPUT_SAM"
done

# Mensagem de conclusão
echo "Conversão concluída. Arquivo '$OUTPUT_SAM' gerado."
