#!/bin/bash

# script: coverage_analysis.sh
# autores: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 22 AGO 2025
# última atualização: 22 AGO 2025
# Uso: Estima a cobertura de um alinhamento usando samtools

#
# Validação da entrada de dados na linha de comando
#
SAMFILENAME=$1	# Caminho/Nome do arquivo.sam a ser analisado

#
# Análise propriamente dita
#

# Converte um arquivo .sam em .bam
samtools view -b -o "${SAMFILENAME}.bam" ${SAMFILENAME}

# Ordena as sequencias antes da análise de cobertura propriamente dita
samtools sort "${SAMFILENAME}.bam" -o "${SAMFILENAME}.sorted.bam"

# Estima a profundidade em cada posição do arquivo .bam
samtools depth ${SAMFILENAME}.sorted.bam" -o ${SAMFILENAME}.sorted.bam.depth"

# Estima a cobertura média do alinhamento
samtools coverage -A -w -r ${SAMFILENAME}.sorted.bam"


