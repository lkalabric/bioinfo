# Script bioinfo_isis.sh
# Autor: Isis Katarina
# Data da criação: 30/09/2026
# Sintáse: bioinfo_isis.sh <AMOSTRA>

# Requisitos:
# - Linux: fastqc, fastp, spades
# - Miniconda: https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install
# curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash ~/Miniconda3-latest-Linux-x86_64.sh
# source ~/.bashrc
# conda create -n pacote_env
# conda activate pacote_env
# conda install pacote
# - Pacotes Conda: trimmomatic, multiqc, quast

# 1) Configuração da entrada de dados
AMOSTRA=$1
# AMOSTRA="102390"
INPUT_DIR="data/hermes"

# Valida se a variável $AMOSTRA está vazia ou se é existente
if [ -z "${AMOSTRA}" ]; then
    echo "Erro: Insira o nome da amostra. Sintáxe: bioinfo-isis.sh 102390"
    echo "Encerrando o script..."
    exit 1
else
    if [ ! $(find "${INPUT_DIR}" -name "${AMOSTRA}*R1*" -print -quit) ]; then
        echo "Erro: Amostra não encontrada! Entrar com um nome de amostra válido."
        echo "Encerrando o script..."
        exit 2
    fi
fi

# 2) Configuração da saída de dados
OUTPUT_DIR="qc-results/${AMOSTRA}"
R1=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name "${AMOSTRA}*R1*" -printf "%f\n" -quit)
R2=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name "${AMOSTRA}*R2*" -printf "%f\n" -quit)
REFSEQ="data/refseq/NC_045512_sequence.fasta" # Referencia de Sars-Cov2
THREADS=$(nproc)

# Caminhos dos diretórios das análises parciais
QC_RAW_DIR="${OUTPUT_DIR}/qc_raw"
QC_RAW_MULTIQC_DIR="${OUTPUT_DIR}/qc_raw_multiqc"
FASTP_DIR="${OUTPUT_DIR}/fastp_out"
TRIM_DIR="${OUTPUT_DIR}/trim_out"
SPADES_DIR="${OUTPUT_DIR}/spades_out"
QUAST_DIR="${OUTPUT_DIR}/quast_out"
ASSEMBLY_DIR="${OUTPUT_DIR}/assembly_out"

# Cria as pastas com as análises parciais
mkdir -p temp "${QC_RAW_DIR}" "${QC_RAW_MULTIQC_DIR}" "${FASTP_DIR}" "${TRIM_DIR}" "${SPADES_DIR}" "${QUAST_DIR}" "${ASSEMBLY_DIR}"

# 3) Controle de qualidade pelo fastqc
fastqc -t "${THREADS}" -o "${QC_RAW_DIR}" "${INPUT_DIR}/${R1}" "${INPUT_DIR}/${R2}"

# 3.1) Controle de qualidade pelo multiqc
# Previne a mensagem "executar conda init primeiro"
# Encontre e carregue a função 'conda' para o subshell
# (Ajuste o caminho para a sua instalação do Conda ou Miniconda)
source ~/miniconda3/etc/profile.d/conda.sh
# Se usar anaconda: source ~/anaconda3/etc/profile.d/conda.sh

# Ativa o ambiente Conda contendo MultiQC
conda activate multiqc
multiqc \
    "${QC_RAW_DIR}" \
    -o "${QC_RAW_MULTIQC_DIR}"

# 4) Pré-processamento dos dados pelo fastp
fastp \
  --in1 "${INPUT_DIR}/${R1}" \
  --in2 "${INPUT_DIR}/${R2}" \
  --out1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" \
  --out2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 30 \
  --length_required 50 \
  --correction \
  --html "${FASTP_DIR}/${AMOSTRA}_fastp_report.html" \
  --json "${FASTP_DIR}/${AMOSTRA}_fastp_report.json" \
  --thread "${THREADS}"

# 4.1) Pré-processamento dos dados pelo Trimmomatic
# Ativa o ambiente Conda contendo o Trimmomatic
conda activate trimmomatic
# Define o caminho do adaptador dinamicamente
ADAPTERS="${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa"
trimmomatic PE -threads "${THREADS}" \
  "${INPUT_DIR}/${R1}" "${INPUT_DIR}/${R2}" \
  "${TRIM_DIR}/${AMOSTRA}_R1_paired.fq.gz" "${TRIM_DIR}/${AMOSTRA}_R1_unpaired.fq.gz" \
  "${TRIM_DIR}/${AMOSTRA}_R2_paired.fq.gz" "${TRIM_DIR}/${AMOSTRA}_R2_unpaired.fq.gz" \
  ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# 5) Montagem de novo usando Spades
spades.py \
    -1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" \
    -2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" \
    -o "${SPADES_DIR}/1" \
    --threads "${THREADS}"

# 5.1) Avaliação da montagem
# Link: https://github.com/ablab/quast
# Link: https://anaconda.org/channels/bioconda/packages/quast/overview
# Instala Quast num ambiente Conda com Python 3.10 ou 3.11 (compatível com a biblioteca padrão distutils)
# conda install quast python=3.10 -y
# Ativa o ambiente Conda contendo o Quast
conda activate quast
# quast.py "${SPADES_DIR}/1/contigs.fasta" -o ${QUAST_DIR}
quast.py \
    "${SPADES_DIR}/1/contigs.fasta" \
    -o "${QUAST_DIR}/1" \
    -r "${REFSEQ}"

# 6) Montagem por referência usando Spades
# Análise usando o preset --metaviral do Spades
spades.py \
    --metaviral \
    -1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" \
    -2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" \
    -o "${SPADES_DIR}/2" \
    --threads "${THREADS}"
# Análise usando uma sequência de referência
spades.py \
    -1 "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" \
    -2 "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" \
    -o "${SPADES_DIR}/3" \
    --trusted-contigs "${REFSEQ}" \
    --threads "${THREADS}"

# 6.1) Avaliação da montagem
# Ativa o ambiente Conda contendo o Quast
conda activate quast
quast.py \
    "${SPADES_DIR}/2/contigs.fasta" \
    -o "${QUAST_DIR}/2" \
    -r "${REFSEQ}"
quast.py \
    "${SPADES_DIR}/3/contigs.fasta" \
    -o "${QUAST_DIR}/3" \
    -r "${REFSEQ}"

# 6.2) Montagem por referência/Mapeamento usando o bwa-mem2
# Link: https://github.com/bwa-mem2/bwa-mem2
# bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam
# Ativa o ambiente Conda contendo o bwa-mem2
conda activate bwa-mem2
bwa-mem2 index \
    "${REFSEQ}"
bwa-mem2 mem \
    -t "${THREADS}" \
    "${REFSEQ}" \
    "${FASTP_DIR}/${AMOSTRA}_R1.clean.fastq.gz" \
    "${FASTP_DIR}/${AMOSTRA}_R2.clean.fastq.gz" \
    | samtools sort -o "${ASSEMBLY_DIR}/alinhado.bam"

exit 3

### Em desenvolvimento

# Relatórios do Mapeamento
samtools flagstat alinhado.bam > relatorio_mapeamento.txt
samtools stats alinhado.bam > estatisticas_detalhadas.txt
samtools coverage alinhado.bam

# Sequência consenso
THREADS=$(nproc)
REF="data/refseq/NC_045512_sequence.fasta"

# 1. Mapear genótipos e variantes (gera um VCF comprimido)
bcftools mpileup -Ou -f "$REF" alinhado.bam | \
bcftools call -mv -Oz -o variantes.vcf.gz

# 2. Indexar o arquivo VCF
bcftools index variantes.vcf.gz

# 3. Gerar a sequência consenso (FASTA final)
bcftools consensus -f "$REF" variantes.vcf.gz > sequencia_consenso.fasta
