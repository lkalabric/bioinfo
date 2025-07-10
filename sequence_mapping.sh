#!/bin/bash

# Mapea sequencias query multiseqfasta file (.fasta) em uma refseq (banco de dados Blast ou Diamond)
# Autor: Luciano Kalabric Silva
# Data de criação: 10/07/2025
# Última atualização: 10/07/2025
# Log: Debugging
# 15/04/2024 - Permite que o script continue de onde parou para permitir concluir a criação do arquivo .map
# Sintáxe: Em revisão

#
# Validação da entrada de dados na linha de comando
#
DBTOOL=$1  # Ferramenta de banco de dados Blast ou Diamond
TAXON=$2   # Taxon path/filename or taxondir
DBNAME=$3  # Database name
DBTYPE=$4  # Tipo de sequencia a ser analisada nucl ou prot
QUERY=$5   # Arquivo .fasta com as sequencias query
if [[ $# -lt 4 ]]; then
	echo "Falta algum parâmetro: (1) ferramenta de geração de banco de dados, (2) caminho/nome do Taxon, (3) apelido do banco de dados a ser criado, (4) tipo do sequencia, ou (5) o nome do arquivo Query (.fasta)!"
	echo "Sintáxe: ./sequence_mapping.sh <-blast/-diamond> <TAXONDIR/TAXONFILENAME> <DBNAME> <BDTYPE: nucl/prot> <QUERY>"
 	exit 0
fi

# Criação do banco de dados
echo "Criando o banco de dados $DBNAME..."
fasta2db.sh $DBTOOL $TAXON $DBNAME $DBTYPE
