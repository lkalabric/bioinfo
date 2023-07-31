#!/bin/bash

# script: blastanything.sh
# autores: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 28 MAR 2023
# última atualização: 28 MAR 2023
# Uso: Faz busca em banco de dados Blast

#
# Validação da entrada de dados na linha de comando
#
QUERY=$1	# Nome do diretório que conte o(s) arquivo(s) query no formato fasta
BLASTDBDIR=$2	# Diretório contendo o BlastDB
BLASTSUITE=$3   # Tipo de busca Blast

# Blast suites disponíveis: 
# blastn - search a nucleotide db using a nucleotide query
# blastp - search a protein db using a protein query
# blastx - search a protein db using a translated nucleotide query
# tblastn - search a translated nucleotide db using a protein query
# tblastx - search a translated nucleotide db using a translated nucleotide query

if [[ $# -lt 3 ]]; then
	echo "Falta o nome do arquivo ou caminho contendo as queries, diretório BlastDB ou blast suite!"
	echo "Sintáxe: ./blastanything.sh <QUERYFILENAME/QUERYDIR> <BLASTDBDIR> <BLASTSUITE: blastn/blastp/blastx/tblastn/tblastx>"
	exit 0
fi

# Saída de dados
BLASTRESULTSDIR="${HOME}/blast-results"
[[ ! -d $BLASTRESULTSDIR ]] && mkdir $BLASTRESULTSDIR

# Preparação do BLASTDB local
# Script: fasta2blastdb.sh
# Concatena todas as REFSEQs num arquivo refseq.fasta único e cria o BLASTDB
# Extrai o label das sequencias em refseq.fasta e cria o arquivo refseq.acc 
# A partir do arquivo refseq.acc, cria o arquivo refseq.map que mapeia os taxid (números que identificam as espécies taxonômica)

# Declaração das funções do script
function blastanything () {
  # Classificação taxonômica das reads utilizando blastn
  # Busca as QUERIES e salva na pasta diretório BLASTNREADSDIR
  # Parâmetros do Blast suite
  # Filter
  EVALUE=0.000001
  QCOV=90
  # outfmt
  # Link: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
  echo -e "Classificando as reads pelo ${BLASTSUITE}...\n"
  if [ -f ${QUERY} ]; then
	# Cria o comando Blast suite para busca em banco de sequencias local
	QUERYNAME=$(basename ${QUERY})
      	CALL_FUNC=echo $(${BLASTSUITE} -db "${BLASTDBDIR}refseq" -query "${QUERY}" -out "${BLASTRESULTSDIR}/${QUERYNAME}.${BLASTSUITE}" -outfmt "6 qseqid sseqid length mismatch gapopen qstart qend sstart send evalue qcovhsp" -qcov_hsp_perc ${QCOV} -max_target_seqs 1)
      	# Executa o comando contido na variável CALL_FUNC
      	eval $CALL_FUNC 
      	# Gera o arquivo de log
	echo "${i} $(wc -l < ${BLASTRESULTSDIR}/${QUERYNAME}.${BLASTSUITE})" >> ${BLASTRESULTSDIR}/passed_reads.log
  else
	  for i in $(find ${QUERY}/*.fasta -type f -exec basename {} .fasta \; | sort); do
		# Cria o comando Blast suite para busca em banco de sequencias local
		CALL_FUNC=echo$(${BLASTSUITE} -db "${BLASTDBDIR}refseq" -query "${QUERY}${i}.fasta" -out "${BLASTRESULTSDIR}/${i}.${BLASTSUITE}" -outfmt "6 qseqid length mismatch gapopen qstart qend sstart send evalue qcovhsp" -qcov_hsp_perc ${QCOV} -max_target_seqs 1)
		# Executa o comando contido na variável CALL_FUNC
		eval $CALL_FUNC 
		# Gera o arquivo de log
		echo "${i} $(wc -l < ${BLASTRESULTSDIR}/${i}.${BLASTSUITE})" >> ${BLASTRESULTSDIR}/passed_reads.log
	  done
  fi
  echo "Resultados ${BLASTSUIE} obtidos com sucesso!"
}


#
# Main do script
#
# Executa a busca usando o Blast suite informado
blastanything
exit 1
