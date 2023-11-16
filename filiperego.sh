
USERNAME=kalabric

cp -r filiperego/qc/ ${USERNAME}/qc/

QUALITY=30
LENTGH=50
HEAD=14
TAIL_R1=1
TAIL_R2=1

cd ${USERNAME}/qc
conda activate fastp
fastp -i 292879835_S26_L001_R1_001.fastq.gz -I 292879835_S26_L001_R2_001.fastq.gz -o 292879835_R1_trimmed.fastq.gz -O 292879835_R2_trimmed.fastq.gz -q ${QUALITY} -l ${LENGTH} -f ${HEAD} -t ${TAIL_R1} -T ${TAIL_R2} -h 292879935.html
conda deactivate

bwa index sars_cov_2.fasta
bwa mem sars_cov_2_ref.fasta 292879835_R1_trimmed.gz 292879835_R2_trimmed.gz | gzip -3 > aln-pe_292879835
