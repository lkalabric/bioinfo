#!/bin/bash

# Cria uma lista de Taxid a partir do número de uma lista de Acc para análise pelo BLAST
# Autor: Luciano Kalabric Silva
# Data: 18/12/2023
# Última atualização: 18/12/2023
# Referencia: https://www.ncbi.nlm.nih.gov/books/NBK569841/

# Sintáxe: 
#	./acc2taxid.sh <TAXONDIR/TAXONFILENAME>"

#
# Validação da entrada de dados na linha de comando
#
TAXON=$1	# Taxon path/filename or taxondir
if [[ $# -lt 1 ]]; then
	echo "Falta o caminho/nome ou o caminho do Taxon!"
	echo "Sintáxe: ./acc2taxid.sh <TAXONDIR/TAXONFILENAME>"
	exit 0
fi

# Diretórios dos dados (Subjects)
# Salvar os arquivos contendo as sequencias referências no formato Fasta (obtidos do Genbank) 
# individualmente ou no formato multiseq Fasta neste diretório ou em sub-diretórios. Os arquivos
# serão contatenados recursivamente em um único arquivo refseq.fasta para criação do banco de dados
ACC2TAXIDDIR="${HOME}/data/TAXID"
mkdir -p "${ACC2TAXIDDIR}"

# Se TAXON for um diretório, concatena todos os arquivos .fasta em ${BLASTDBNAME}.fasta
# Exceto o arquivo ${BLASTDBNAME}/refseq.fasta que é gerado pelo make_refgen.sh
echo "Concatenando as sequencias referências em ${TAXON}..."
if [ -f ${TAXON} ]; then
	cat ${TAXON} > "${ACC2TAXIDDIR}/taxons.fasta"
else
	# find ${TAXON} -type f -iname '*.fasta' -print0 | sort -z | xargs -0 cat > "${BLASTDBDIR}/refseq.fasta"
	find ${TAXON} -name '*.fasta' -exec cat {} + > "${ACC2TAXIDDIR}/taxons.fasta"
fi

# Processa a linha de descrição das sequencias referências para conter apenas o número de acesso sem espaços
echo "Processando os labels do arquivo ${ACC2TAXIDDIR}/taxons.fasta..."
[[ -f ${ACC2TAXIDDIR}/taxons.old ]] && rm ${ACC2TAXIDDIR}/taxons.old
mv ${ACC2TAXIDDIR}/taxons.fasta ${ACC2TAXIDDIR}/taxons.old
# Cria a lista de números de acc Genbank a partir do arquivo .fasta
echo "Criando o arquivo taxons.acc..."
grep ">" ${ACC2TAXIDDIR}/taxons.old | cut -d "." -f 1 | cut -c 2-10 > ${ACC2TAXIDDIR}/taxons.acc

# Cria a lista de taxid a partir nos números de acc Genbank
[[ -f ${ACC2TAXIDDIR}/taxons.map ]] && rm ${ACC2TAXIDDIR}/taxons.map
# Retrive Taxid
echo "Criando o arquivo ${ACC2TAXIDDIR}/taxons.map..."
while read -r line; do
  # echo "$line "$(efetch -db nuccore -id "$line" -format docsum | xtract -pattern DocumentSummary -element TaxId) >>${BLASTDBDIR}/refseq.map
  # Alternativamente, podemos obter o Taxid usado esearch em combinação com esummary
  echo "Buscando acc $line"
  esearch -db assembly -q $line | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Taxid >> ${ACC2TAXIDDIR}/taxons.map
	# echo "$line $(esearch -db assembly -q "$line" | esummary | xtract -pattern DocumentSummary -element Taxid)" >> ${ACC2TAXIDDIR}/taxons.map
done < ${ACC2TAXIDDIR}/taxons.acc
