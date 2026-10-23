# Pacote requerido: bedtools
# sudo apt install bedtools

# Arquivo .bed a ser convertido pelo script
ARQUIVO=$1

if [[ -z "${ARQUIVO}" || ! -f "${ARQUIVO}" ]]; then
    echo "Uso: bed2fasta.sh <ARQUIVO>."
    echo "Arquivo ausente ou não encontrado. Encerrando o script..."
    exit 1
fi


# Sequência de referência utilizada pelo kit Illumina
REFSEQ="MN908947.3.fasta"

# Análise do arquivo .bed
awk '{print NF}' ${ARQUIVO} | sort | uniq -c

# Separa os primers forward e reverse em dois .bed
grep "LEFT" ${ARQUIVO} > primer_forward.bed
grep "RIGHT" ${ARQUIVO} > primer_reverse.bed

# Converte .bed para .fasta usando o comando bedtools getfasta
bedtools getfasta -fi ${REFSEQ} -bed primer_forward.bed -s -name -fo primer_forward.fasta
bedtools getfasta -fi ${REFSEQ} -bed primer_reverse.bed -s -name -fo primer_reverse.fasta
