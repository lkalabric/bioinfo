#!/bin/bash

# Autor: Filipe Rego
# Uso: Quality control

USERNAME=kalabric
#scp -r bioinfo@200.128.7.33:filiperego/qc/ ${USERNAME}/qc/
#scp -r bioinfo@192.168.34.4:filiperego/qc/ ${USERNAME}/qc/
cp -r filiperego/qc/ ${USERNAME}/qc/

FILENAME="292879835_S26_L001"
QUALITY=30
LENGTH=50
HEAD=14
TAIL_R1=1
TAIL_R2=1

cd ${USERNAME}/qc
source activate fastp
fastp -i ${FILENAME}_R1_001.fastq.gz -I ${FILENAME}_R2_001.fastq.gz -o ${FILENAME}_R1_trimmed.fastq.gz -O ${FILENAME}_R2_trimmed.fastq.gz -q ${QUALITY} -l ${LENGTH} -f ${HEAD} -t ${TAIL_R1} -T ${TAIL_R2} -h 292879935.html
conda deactivate

# Cria o arquivo de índice para mapeamento da refseq
bwa index sars_cov_2_ref.fasta

# Monta por referência das reads
bwa mem sars_cov_2_ref.fasta ${FILENAME}_R1_trimmed.fastq.gz ${FILENAME}_R2_trimmed.fastq.gz | gzip -3 > aln-pe_${FILENAME}

exit 0

#gere a sequencia consenso
#transforma o bam em sorted bam
samtools view -bS aln-pe_${FILENAME} | samtools sort - -o ${FILENAME}.bam

#cria a consenso fastq
samtools mpileup -uf sars_cov_2_ref.fasta ${FILENAME}.bam | bcftools call -c | vcfutils.pl vcf2fq > ${FILENAME}.fastq

#cria a consenso fasta (não faremos pois esqueci de pedir para instalar o programa)
seqtk seq -aQ64 -q20 -n ${FILENAME}.fastq > ${FILENAME}.fasta

#importe a sequencia consenso no formato fastq e abra no aliview
