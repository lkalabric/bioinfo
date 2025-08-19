#!/bin/bash

# Este script realiza o mapeamento de sequências de leituras curtas (FASTQ)
# em um genoma de referência (FASTA) usando BWA e SAMtools.

# Uso:
# ./seq_mapping.sh <genoma_referencia.fasta> <leituras_R1.fastq> <leituras_R2.fastq>

#---------------------------------------------------------------------------------------------------
# 1. VALIDAÇÃO E DEFINIÇÃO DE VARIÁVEIS
#---------------------------------------------------------------------------------------------------

# Certifica que o script irá parar se qualquer comando falhar
set -e

# Verifica se os três argumentos obrigatórios foram fornecidos
if [ "$#" -lt 3 ]; then
  echo "Uso: $0 <genoma_referencia.fasta> <leituras_R1.fastq> <leituras_R2.fastq>"
  echo "Exemplo: $0 Ecoli_K12.fasta reads_1.fastq.gz reads_2.fastq.gz"
  exit 1
fi

# Define variáveis com base nos argumentos da linha de comando
referencia="$1"
leituras_R1="$2"
leituras_R2="$3"

# Define o nome dos arquivos de saída com base no nome do arquivo de referência
nome_base_ref=$(basename "$referencia" .fasta)
saida_bam_sorted="${nome_base_ref}_mapeado_ordenado.bam"
saida_log="${nome_base_ref}_mapeamento.log"

echo "Iniciando o pipeline de mapeamento..."
echo "-------------------------------------"
echo "Genoma de referência: $referencia"
echo "Leituras R1: $leituras_R1"
echo "Leituras R2: $leituras_R2"
echo "Arquivo de saída: $saida_bam_sorted"
echo "Log: $saida_log"
echo "-------------------------------------"

#---------------------------------------------------------------------------------------------------
# 2. INDEXAÇÃO DA REFERÊNCIA (FASTA)
#    Esta etapa só precisa ser feita uma vez para cada genoma de referência.
#---------------------------------------------------------------------------------------------------

if [ ! -f "$referencia.fai" ]; then
  echo "Iniciando a indexação do genoma de referência..." | tee -a "$saida_log"
  bwa index "$referencia" >> "$saida_log" 2>&1
  echo "Indexação concluída." | tee -a "$saida_log"
else
  echo "Genoma de referência já indexado. Pulando a indexação." | tee -a "$saida_log"
fi

#---------------------------------------------------------------------------------------------------
# 3. MAPEAMENTO (ALIGNMENT) USANDO BWA E PÓS-PROCESSAMENTO COM SAMTOOLS
#---------------------------------------------------------------------------------------------------

echo "Iniciando o mapeamento das leituras com BWA..." | tee -a "$saida_log"

bwa mem "$referencia" "$leituras_R1" "$leituras_R2" 2>> "$saida_log" | \
samtools view -bS - 2>> "$saida_log" | \
samtools sort -o "$saida_bam_sorted" - 2>> "$saida_log"

echo "Mapeamento e pós-processamento concluídos." | tee -a "$saida_log"

#---------------------------------------------------------------------------------------------------
# 4. INDEXAÇÃO DO ARQUIVO BAM
#    Necessário para visualização em genomas browsers (como o IGV).
#---------------------------------------------------------------------------------------------------

echo "Iniciando a indexação do arquivo BAM..." | tee -a "$saida_log"
samtools index "$saida_bam_sorted" >> "$saida_log" 2>&1
echo "Indexação do BAM concluída." | tee -a "$saida_log"

echo "-------------------------------------"
echo "Pipeline finalizado com sucesso."
echo "Arquivos de saída:"
echo "-> Alinhamento ordenado e indexado: $saida_bam_sorted"
echo "-> Índice do BAM: ${saida_bam_sorted}.bai"
echo "-> Log do processo: $saida_log"
