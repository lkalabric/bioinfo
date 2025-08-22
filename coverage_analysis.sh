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

# Convete uma aquivo .sam em .bam
samtools view -b -o output.bam output.sam

# Estima a profundidade em cada posição do arquivo .bam
samtools depth output.bam -o output.bam.depth

