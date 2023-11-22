#!/bin/bash

# Autor: Luciano Kalabric
# Uso: BLAST CLI

cd

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Número de parametros ilegal."
    echo "Síntaxe: kalabric_blast.sh <username_dir> <dataset_contig.fasta>"
    echo "Exemplo: kalabric_blast.sh kalabric dataset1-contig1.fasta"
    exit 0    
fi
USERNAME_DIR=$1
CONTIG=$2

if [ ! -d "${HOME}/${USERNAME_DIR}" ]; then
    echo "Directório ${USERNAME_DIR}/ não existe."
    exit 0
else
    if [ ! -d "${HOME}/${USERNAME_DIR}/blast-analysis" ]; then
        # Cria a árvore de diretórios e copia os dados para realização das análises
        echo "Criando árvore de diretórios para análise BLAST CLI..."
        cd ${USERNAME_DIR}
        mkdir blast-analysis
        cd blast-analysis/
        mkdir blastdb
        mkdir queries
        mkdir results
        # Copia os dados de exemplo para a pasta queries/
        cp ${HOME}/examples/blast/*.fasta queries/
    else
        echo "Árvore de diretórios criada com sucesso..."
        cd ${USERNAME_DIR}
        cd blast-analysis/
        cp ${HOME}/examples/blast/*.fasta queries/
    fi
fi

# Cria o banco de dados apenas uma vez
if [ ! -f "${HOME}/${USERNAME_DIR}/blast-analysis/blastdb/virome.nhr" ]; then
    echo "Criando banco de dados virome..."
    makeblastdb -in ~/examples/blast/viral.1.1.genomic.fna -dbtype nucl -out "${HOME}/${USERNAME_DIR}/blast-analysis/blastdb/virome"
    ls blastdb
else
    echo "Banco de dados virome já criado!"
fi

# Realiza a busca por similaridade utilizando BLASTN e retorna report do apenas do best hit (e-value <= 1E-6, cobertura >= 90%, formato tabular)
echo "Realizando a busca da query ${CONTIG} no banco de sequências virome..." 
blastn -db blastdb/virome -query queries/${CONTIG} -out results/${CONTIG}.e-6c90hsp1.blastn -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1 -outfmt "6 sacc staxid"
cat results/${CONTIG}.e-6c90hsp1.blastn

