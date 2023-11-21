#!/bin/bash

# Autor: Luciano Kalabric
# Uso: BLAST CLI

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Número de parametros ilegal."
    echo "Síntaxe: kalabric_blast.sh <username_dir> <dataset_contig.fasta>"
    echo "Exemplo: kalabric_blast.sh kalabric dataset1_contig1.fasta"
    exit 0    
fi
USERNAME_DIR=$1
CONTIG=$2

if [ ! -d "${HOME}/${USERNAME_DIR}" ]; then
    echo "Directório ${USERNAME_DIR}/ não existe."
    exit 0
else
    # Cria a árvore de diretórios e copia os dados para realização das análises
    cd "${USERNAME_DIR}"
    mkdir blast-analysis
    cd blast-analysis/
    mkdir blastdb
    mkdir queries
    mkdir results
    # Copia os dados de exemplo para a pasta queries/
    cp "${HOME}/examples/blast/*.fasta" "${USERNAME_DIR}/blast-analysis/queries"
fi

# Cria o banco de dados
echo "Criando banco de dados virome..."
makeblastdb -in ~/examples/blast/viral.1.1.genomic.fna -dbtype nucl -out "${HOME}/${USERNAME_DIR}/blast-analysis/blastdb/virome"
ls blastdb

# Realiza a busca por similaridade utilizando BLASTN e retorna report do apenas do best hit (e-value <= 1E-6, cobertura >= 90%, formato tabular)
blastn -db blastdb/viral.1.1.genomic -query queries/${CONTIG} -out results/${CONTIG}.blastn -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1 -outfmt "6 sacc staxid"
