#!/bin/bash

# Autor: Filipe Rego & Luciano Kalabric
# Uso: Quality control

cd

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Número de parametros ilegal."
    echo "Síntaxe: filiperego_qc.sh <username_dir> <runname>"
    echo "Examplo: filiperego_qc.sh kalabric 292879835_S26_L001"
    exit 0    
fi
USERNAME_DIR=$1
RUNNAME=$2

if [ ! -d "${HOME}/${USERNAME_DIR}" ]; then
    echo "Directório ${USERNAME_DIR}/ não existe."
    exit 0
else
    cd ${USERNAME_DIR}
    mkdir qc/
    cp ${HOME}/examples/ngs-qc/* qc/
    cd qc/
fi

#RUNNAME="292879835_S26_L001"
QUALITY=30
LENGTH=50
HEAD=14
TAIL_R1=1
TAIL_R2=1

source activate fastp
fastp -i ${RUNNAME}_R1_001.fastq.gz -I ${RUNNAME}_R2_001.fastq.gz -o ${RUNNAME}_R1_trimmed.fastq.gz -O ${RUNNAME}_R2_trimmed.fastq.gz -q ${QUALITY} -l ${LENGTH} -f ${HEAD} -t ${TAIL_R1} -T ${TAIL_R2} -h 292879935.html
conda deactivate

# Cria o arquivo de índice para mapeamento da refseq
bwa index sars_cov_2_ref.fasta

# Montagem por referência das reads
bwa mem sars_cov_2_ref.fasta ${RUNNAME}_R1_trimmed.fastq.gz ${RUNNAME}_R2_trimmed.fastq.gz | gzip -3 > aln-pe_${RUNNAME}.sam.gz

# Gere a sequencia consenso
#transforma o sam->bam em sorted bam
samtools view -bS aln-pe_${RUNNAME}.sam.gz | samtools sort - -o ${RUNNAME}.bam
# cria um indice para o arquivo bam
samtools index ${RUNNAME}.bam ${RUNNAME}.bam.bai

#cria a consenso fastq
samtools mpileup -uf sars_cov_2_ref.fasta ${RUNNAME}.bam | bcftools call -c | vcfutils.pl vcf2fq > ${RUNNAME}.fastq

#cria a consenso fasta (não faremos pois esqueci de pedir para instalar o programa)
seqtk seq -aQ64 -q20 -n ${RUNNAME}.fastq > ${RUNNAME}.fasta

#importe a sequencia consenso no formato fastq e abra no aliview
