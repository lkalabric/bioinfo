# Script bioingo_isis.sh
# Autor: Isis Katarina
# Data da criação: 30/09/2026
# Sintáse: bioinfo_isis.sh <AMOSTRA>

AMOSTRA=$1
# Valida se a variável está vazia
if [ -z "$AMOSTRA" ]; then
    echo "Erro: Amostra inexistente ou nome vazio. Entrar com um nome de amostra válido."
    echo "Encerrando o script..."
    exit 1
fi

# Salvar os resultados das análises por amostra no diretório results/
INPUT_DIR="data/hermes"
# AMOSTRA="102390"
OUTPUT_DIR="qc-results/${AMOSTRA}"
R1=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "${AMOSTRA}*R1*" -printf "%f\n" -quit)
R2=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "${AMOSTRA}*R2*" -printf "%f\n" -quit)

# Caminhos dos diretórios das análises parciais
QC_RAW_DIR="${OUTPUT_DIR}/qc_raw"
QC_RAW_MULTIQC_DIR="${OUTPUT_DIR}/qc_raw_multiqc"
FASTP_DIR="${OUTPUT_DIR}/fastp_out"
TRIM_DIR="${OUTPUT_DIR}/trim_out"

# Cria as pastas com as análises parciais
mkdir -p temp "$QC_RAW_DIR" "$QC_RAW_MULTIQC_DIR" "$FASTP_DIR" "$TRIM_DIR"

# Controle de qualidade
fastqc -t 4 -o "$QC_RAW_DIR" "$INPUT_DIR/$R1" "$INPUT_DIR/$R2"

# Previne a mensagem "executar conda init primeiro"
# Encontre e carregue a função 'conda' para o subshell
# (Ajuste o caminho para a sua instalação do conda ou miniconda)
source ~/miniconda3/etc/profile.d/conda.sh
# Se usar anaconda: source ~/anaconda3/etc/profile.d/conda.sh

# Ativa o ambiente normalmente
conda activate multiqc
multiqc "$QC_RAW_DIR" -o "$QC_RAW_MULTIQC_DIR"

# Pré-processamento dos dados
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

trimmomatic PE -threads 4 \
  "$TRIM_DIR/$R1" "$TRIM_DIR/$R2" \
  "${TRIM_DIR}/${AMOSTRA}_R1_paired.fq.gz" "${TRIM_DIR}/${AMOSTRA}_R1_unpaired.fq.gz" \
  "${TRIM_DIR}/${AMOSTRA}_R2_paired.fq.gz" "${TRIM_DIR}/${AMOSTRA}_R2_unpaired.fq.gz" \
  ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:True \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
