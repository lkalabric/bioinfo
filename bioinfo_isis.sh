# Script bioinfo_isis.sh
# Autor: Isis Katarina
# Data da criação: 30/09/2026
# Sintáse: bioinfo_isis.sh <AMOSTRA>

# Requisitos:
# - Linux: fastqc, fastp, spades
# - Conda: trimmomatic, multiqc, quast

# Nome da amostra passada na linha de comando
AMOSTRA=$1
# 1) Valida se a variável está vazia
if [ -z "$AMOSTRA" ]; then
    echo "Erro: Amostra inexistente ou nome vazio. Entrar com um nome de amostra válido."
    echo "Encerrando o script..."
    exit 1
fi

# 2) Entrada e saída de dados
INPUT_DIR="data/hermes"
# AMOSTRA="102390"
OUTPUT_DIR="qc-results/${AMOSTRA}"
R1=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "${AMOSTRA}*R1*" -printf "%f\n" -quit)
R2=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "${AMOSTRA}*R2*" -printf "%f\n" -quit)
REFSEQ="data/refseq/NC_045512_sequence.fasta" # Referencia de Sars-Cov2

# Caminhos dos diretórios das análises parciais
QC_RAW_DIR="${OUTPUT_DIR}/qc_raw"
QC_RAW_MULTIQC_DIR="${OUTPUT_DIR}/qc_raw_multiqc"
FASTP_DIR="${OUTPUT_DIR}/fastp_out"
TRIM_DIR="${OUTPUT_DIR}/trim_out"
SPADES_DIR="${OUTPUT_DIR}/spades_out"
QUAST_DIR="${OUTPUT_DIR}/quast_out"
ASSEMBLY_DIR="${OUTPUT_DIR}/assembly_out"

# Cria as pastas com as análises parciais
mkdir -p temp "$QC_RAW_DIR" "$QC_RAW_MULTIQC_DIR" "$FASTP_DIR" "$TRIM_DIR" "$SPADES_DIR" "$QUAST_DIR" "$ASSEMBLY_DIR"

# 3) Controle de qualidade pelo fastqc
fastqc -t 4 -o "$QC_RAW_DIR" "$INPUT_DIR/$R1" "$INPUT_DIR/$R2"

# 3.1) Controle de qualidade pelo multiqc
# Previne a mensagem "executar conda init primeiro"
# Encontre e carregue a função 'conda' para o subshell
# (Ajuste o caminho para a sua instalação do Conda ou Miniconda)
source ~/miniconda3/etc/profile.d/conda.sh
# Se usar anaconda: source ~/anaconda3/etc/profile.d/conda.sh

# Ativa o ambiente Conda contendo MultiQC
conda activate multiqc
multiqc "$QC_RAW_DIR" -o "$QC_RAW_MULTIQC_DIR"

# 4) Pré-processamento dos dados pelo fastp
fastp \
  --in1 "$INPUT_DIR/$R1" \
  --in2 "$INPUT_DIR/$R2" \
  --out1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" \
  --out2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 30 \
  --length_required 50 \
  --correction \
  --html "${FASTP_DIR}/${AMOSTRA}_fastp_report.html" \
  --json "${FASTP_DIR}/${AMOSTRA}_fastp_report.json" \
  --thread 4

# 4.1) Pré-processamento dos dados pelo Trimmomatic
# Ativa o ambiente Conda contendo o Trimmomatic
conda activate trimmomatic
# Define o caminho do adaptador dinamicamente
ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"
trimmomatic PE -threads 4 \
  "$INPUT_DIR/$R1" "$INPUT_DIR/$R2" \
  "${TRIM_DIR}/${AMOSTRA}_R1_paired.fq.gz" "${TRIM_DIR}/${AMOSTRA}_R1_unpaired.fq.gz" \
  "${TRIM_DIR}/${AMOSTRA}_R2_paired.fq.gz" "${TRIM_DIR}/${AMOSTRA}_R2_unpaired.fq.gz" \
  ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# 5) Montagem de novo usando Spades
spades.py -1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" -2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" -o "${SPADES_DIR}/1" --threads 8

# 5.1) Avaliação da montagem
# Link: https://github.com/ablab/quast
# Link: https://anaconda.org/channels/bioconda/packages/quast/overview
# Instala Quast num ambiente Conda com Python 3.10 ou 3.11 (compatível com a biblioteca padrão distutils)
# conda install quast python=3.10 -y
# Ativa o ambiente Conda contendo o Quast
conda activate quast
# quast.py "${SPADES_DIR}/1/contigs.fasta" -o ${QUAST_DIR}
quast.py "${SPADES_DIR}/1/contigs.fasta" -o "${QUAST_DIR}/1" -r "${REFSEQ}"

# 6) Montagem por referência usando Spades
# Análise usando o preset --metaviral do Spades
spades.py --metaviral -1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" -2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" -o "${SPADES_DIR}/2" --threads 8
# Análise usando uma sequência de referência
spades.py -1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" -2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" -o "${SPADES_DIR}/3" --trusted-contigs "${REFSEQ}" --threads 8

# 6.1) Avaliação da montagem
# Ativa o ambiente Conda contendo o Quast
conda activate quast
quast.py "${SPADES_DIR}/2/contigs.fasta" -o "${QUAST_DIR}/2" -r "${REFSEQ}"
quast.py "${SPADES_DIR}/3/contigs.fasta" -o "${QUAST_DIR}/3" -r "${REFSEQ}"

# 6.2) Montagem por referência/Mapeamento usando o bwa-mem2
# Link: https://github.com/bwa-mem2/bwa-mem2
# bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam
# Ativa o ambiente Conda contendo o bwa-mem2
conda activate bwa-mem2
bwa-mem2 mem -t 8 "${REFSEQ}" "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" | samtools sort -o "${ASSEMBLY_DIR}/alinhado.bam"

