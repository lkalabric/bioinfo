#!/bin/bash

# Mapear sequencias query multiseqfasta file (.fasta) em uma refseq
# Autor: Luciano Kalabric Silva
# Data de criação: 19/08/2025
# Última atualização: 19/08/2025
# Log: Debugging
# Sintáxe: Em revisão

#
# Validação da entrada de dados na linha de comando
#
MAPPINGTOOL=$1  # Ferramenta de mapeamento minimap2 ou bwa
REFSEQ=$2   # path/filename refseq.fasta
DBTYPE=$4  # Tipo de sequencia a ser analisada nucl ou prot
QUERY=$5   # path/filename query.fasta com as sequencias query
if [[ $# -lt 4 ]]; then
	echo "Falta algum parâmetro: (1) ferramenta de geração de banco de dados, (2) caminho/nome do Taxon, (3) apelido do banco de dados a ser criado, (4) tipo do sequencia, ou (5) o nome do arquivo Query (.fasta)!"
	echo "Sintáxe: ./sequence_mapping.sh <MAPPINGTOOL: minimpa2/bwa> <REFSEQ_DIR/FILENAME.fasta> <BDTYPE: nucl/prot> <QUERY.fasta>"
 	exit 0
fi

#
# Uso do minimap2
#
# Link: https://github.com/lh3/minimap2

# Instala pacote caso não exista
PACKAGE_NAME="minimap2" # Replace with the actual package name
COMMAND_NAME="minimap2" # Replace with a command provided by the package

if ! which "$COMMAND_NAME" > /dev/null; then
    echo "$PACKAGE_NAME not found. Installing..."
    sudo apt update
    sudo apt install -y "$PACKAGE_NAME"
else
    echo "$PACKAGE_NAME is already installed."
fi

exit


# Baixando e instalando a versão mais recente do minimap2
cd Downloads && git clone https://github.com/lh3/minimap2
cd minimap2 && make
cp minimap2 ~/bin
# or
conda create --name minimap2 python=3.9
conda activate
conda install bioconda::minimap2



# Criação do banco de dados
echo "Criando o banco de dados $DBNAME..."
fasta2db.sh $DBTOOL $REFSEQ $DBTYPE $QUERY



test_dir="examples/minimap2_mapping/"

# Mapeamento propriamente dito
# 1 - long sequences against a reference genome
minimap2 -a $test_dir/MT-human.fa $test_dir/MT-orang.fa > $test_dir/test.sam
exit

# 2 - create an index first and then map
./minimap2 -x map-ont -d MT-human-ont.mmi test/MT-human.fa
./minimap2 -a MT-human-ont.mmi test/MT-orang.fa > test.sam
# 3 - use presets (no test data)
./minimap2 -ax map-pb ref.fa pacbio.fq.gz > aln.sam       # PacBio CLR genomic reads
./minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam         # Oxford Nanopore genomic reads
./minimap2 -ax map-hifi ref.fa pacbio-ccs.fq.gz > aln.sam # PacBio HiFi/CCS genomic reads (v2.19+)
./minimap2 -ax lr:hq ref.fa ont-Q20.fq.gz > aln.sam       # Nanopore Q20 genomic reads (v2.27+)
./minimap2 -ax sr ref.fa read1.fa read2.fa > aln.sam      # short genomic paired-end reads
./minimap2 -ax splice ref.fa rna-reads.fa > aln.sam       # spliced long reads (strand unknown)
./minimap2 -ax splice -uf -k14 ref.fa reads.fa > aln.sam  # noisy Nanopore direct RNA-seq
./minimap2 -ax splice:hq -uf ref.fa query.fa > aln.sam    # PacBio Kinnex/Iso-seq (RNA-seq)
./minimap2 -ax splice --junc-bed=anno.bed12 ref.fa query.fa > aln.sam  # use annotated junctions
./minimap2 -ax splice:sr ref.fa r1.fq r2.fq > aln.sam     # short-read RNA-seq (v2.29+)
./minimap2 -ax splice:sr -j anno.bed12 ref.fa r1.fq r2.fq > aln.sam
./minimap2 -cx asm5 asm1.fa asm2.fa > aln.paf             # intra-species asm-to-asm alignment
./minimap2 -x ava-pb reads.fa reads.fa > overlaps.paf     # PacBio read overlap
./minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf    # Nanopore read overlap
# man page for detailed command line options
man ./minimap2.1

