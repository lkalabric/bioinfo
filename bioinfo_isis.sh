# Salvar os resultados das análises por amostra no diretório results/
INPUT_DIR="data/hermes"
AMOSTRA="102390"
OUTPUT_DIR="qc_results/${AMOSTRA}"
R1="102390_S5_L001_R1_001.fastq.gz"
R2="102390_S5_L001_R2_001.fastq.gz"

# Caminhos dos diretórios das análises parciais
QC_RAW_DIR="${OUTPUT_DIR}/qc_raw"
QC_RAW_MULTIQC_DIR="${OUTPUT_DIR}/qc_raw_multiqc"
FASTP_DIR="${OUTPUT_DIR}/fastp_out"

mkdir temp/              # Cria uma pasta temporária das análises
mkdir "$OUTPUT_DIR"      # Cria a pasta de resultados por amostra
mkdir "$QC_RAW_DIR"          # Cria a pasta para os resultados do FastQC
mkdir "$QC_RAW_MULTIQC_DIR"  # Cria a pasta para os resultados do MultiQC
mkdir "$FASTP_DIR"           # Cria a pasta para os resultados do Fastp

# Controle de qualidade
fastqc -t 4 -o "$QC_RAW" "$INPUT_DIR/$R1" "$INPUT_DIR/$R2"

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
  # --json sample_fastp_report.json \
  --thread 4
