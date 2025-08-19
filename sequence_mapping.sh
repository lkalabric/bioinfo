#!/bin/bash

# Mapear sequencias de nuc em um arquivo query (.fa/.fq) em uma refseq (.fa)
# Autor: Luciano Kalabric Silva
# Data de criação: 19/08/2025
# Última atualização: 19/08/2025
# Log: Debugging
# Sintáxe: sequence_mapping.sh <MAPPINGTOOL: minimap2/bwa/bowtie> <DIR/FILENAME.fasta> <BDTYPE: nucl/prot> <QUERY.fasta ou .fastq>

# Pacotes requeridos: minimap2, bwa
declare -a PACKAGES=("minimap2" "bwa" "bowtie")

# Instala pacote(s) Linux caso não exista(m)
for pack in "${PACKAGES[@]}"; do
	PACKAGE_NAME=$pack # Replace with the actual package name
	COMMAND_NAME=$pack # Replace with a command provided by the package
	if ! which "$COMMAND_NAME" > /dev/null; then
	    echo "$PACKAGE_NAME not found. Installing..."
	    sudo apt update
	    sudo apt install -y "$PACKAGE_NAME"
	else
	    echo "$PACKAGE_NAME is already installed."
	fi
done

#
# Uso do minimap2
#
# Link: https://github.com/lh3/minimap2
# Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database. Typical use cases include: (1) mapping PacBio or Oxford Nanopore genomic reads to the human genome; (2) finding overlaps between long reads with error rate up to ~15%; (3) splice-aware alignment of PacBio Iso-Seq or Nanopore cDNA or Direct RNA reads against a reference genome; (4) aligning Illumina single- or paired-end reads; (5) assembly-to-assembly alignment; (6) full-genome alignment between two closely related species with divergence below ~15%

# Baixando e instalando a versão mais recente do minimap2 diretamente do site do app
# cd Downloads && git clone https://github.com/lh3/minimap2
# cd minimap2 && make
# cp minimap2 ~/bin
# or
# conda create --name minimap2 python=3.9
# conda activate minimap2
# conda install bioconda::minimap2

# Validação da entrada de dados na linha de comando sequence_mapping.sh
MAPPINGTOOL=$1  # Ferramenta de mapeamento minimap2 ou bwa
REFSEQ=$2   # path/filename refseq.fasta
DBTYPE=$3  # Tipo de sequencia a ser analisada nucl ou prot
QUERY=$4   # path/filename query.fasta com as sequencias query
#if [[ $# -lt 4 ]]; then
#	echo "Falta algum parâmetro: (1) ferramenta de geração de banco de dados, (2) caminho/nome do Taxon, (3) apelido do banco de dados a ser criado, (4) tipo do sequencia, ou (5) o nome do arquivo Query (.fasta)!"
#	echo "Sintáxe: ./sequence_mapping.sh <MAPPINGTOOL: minimap2/bwa> <REFSEQ_DIR/FILENAME.fasta> <BDTYPE: nucl/prot> <QUERY.fasta>"
# 	exit 0
#fi

test_dir="examples/minimap2_mapping/"
# Mapeamento propriamente dito
# 1 - long sequences against a reference genome
minimap2 -a $test_dir/MT-human.fa $test_dir/MT-orang.fa > $test_dir/test.sam
# 2 - create an index first and then map
minimap2 -x map-ont -d MT-human-ont.mmi $test_dir/MT-human.fa # indexing minimap2 index format (.mmi)
minimap2 -a MT-human-ont.mmi $test_dir/MT-orang.fa > test.sam # alignment -a generate SAM format output (.sam)
exit
# 3 - use presets (no test data)
# minimap2 -ax map-pb $test_dir/ref.fa $test_dir/pacbio.fq.gz > $test_dir/aln.sam       # PacBio CLR genomic reads
# minimap2 -ax map-ont $test_dir/ref.fa $test_dir/ont.fq.gz > $test_dir/aln.sam         # Oxford Nanopore genomic reads
# minimap2 -ax map-hifi $test_dir/ref.fa $test_dir/pacbio-ccs.fq.gz > $test_dir/aln.sam # PacBio HiFi/CCS genomic reads (v2.19+)
# minimap2 -ax lr:hq $test_dir/ref.fa $test_dir/ont-Q20.fq.gz > $test_dir/aln.sam       # Nanopore Q20 genomic reads (v2.27+)
# minimap2 -ax sr $test_dir/ref.fa $test_dir/read1.fa $test_dir/read2.fa > aln.sam      # short genomic paired-end reads
# minimap2 -ax splice $test_dir/ref.fa $test_dir/rna-reads.fa > $test_dir/aln.sam       # spliced long reads (strand unknown)
# minimap2 -ax splice -uf -k14 $test_dir/ref.fa $test_dir/reads.fa > $test_dir/aln.sam  # noisy Nanopore direct RNA-seq
# minimap2 -ax splice:hq -uf $test_dir/ref.fa $test_dir/query.fa > $test_dir/aln.sam    # PacBio Kinnex/Iso-seq (RNA-seq)
# minimap2 -ax splice --junc-bed=anno.bed12 $test_dir/ref.fa $test_dir/query.fa > $test_dir/aln.sam  # use annotated junctions
# minimap2 -ax splice:sr $test_dir/ref.fa $test_dir/r1.fq $test_dir/r2.fq > $test_dir/aln.sam     # short-read RNA-seq (v2.29+)
# minimap2 -ax splice:sr -j anno.bed12 $test_dir/ref.fa $test_dir/r1.fq $test_dir/r2.fq > $test_dir/aln.sam
# minimap2 -cx asm5 $test_dir/asm1.fa $test_dir/asm2.fa > $test_dir/aln.paf             # intra-species asm-to-asm alignment
# minimap2 -x ava-pb $test_dir/reads.fa $test_dir/reads.fa > $test_dir/overlaps.paf     # PacBio read overlap
# minimap2 -x ava-ont $test_dir/reads.fa $test_dir/reads.fa > $test_dir/overlaps.paf    # Nanopore read overlap
# man page for detailed command line options
man ./minimap2.1

#
# Uso do bwa
#
# Link: https://github.com/lh3/bwa
# BWA is a software package for mapping DNA sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to a few megabases. BWA-MEM and BWA-SW share similar features such as the support of long reads and chimeric alignment, but BWA-MEM, which is the latest, is generally recommended as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

#
# Uso do bowtie
#
# Link: https://github.com/BenLangmead/bowtie2
# Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the genome with an FM Index to keep its memory footprint small: for the human genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment modes.

