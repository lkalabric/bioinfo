
USERNAME=kalabric

cp -r filiperego/qc/ ${USERNAME}/qc/

FILENAME="292879835_S26_L001_"
QUALITY=30
LENTGH=50
HEAD=14
TAIL_R1=1
TAIL_R2=1

cd ${USERNAME}/qc
conda activate fastp
fastp -i ${FILENAME}R1_001.fastq.gz -I ${FILENAME}R2_001.fastq.gz -o ${FILENAME}R1_trimmed.fastq.gz -O ${FILENAME}R2_trimmed.fastq.gz -q ${QUALITY} -l ${LENGTH} -f ${HEAD} -t ${TAIL_R1} -T ${TAIL_R2} -h 292879935.html
conda deactivate

bwa index sars_cov_2_ref.fasta
bwa mem sars_cov_2_ref.fasta ${FILENAME}R1_trimmed.gz ${FILENAME}R2_trimmed.gz | gzip -3 > aln-pe_292879835
