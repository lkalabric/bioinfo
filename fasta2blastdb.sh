#!/bin/bash

# Cria o banco de dados BLAST a partir de arquivos .fasta para busca
# Autor: Luciano Kalabric Silva
# Data: 24/03/2023
# Última atualização: 27/03/2023
# Referencia: https://www.ncbi.nlm.nih.gov/books/NBK569841/
# Log: Debugging

# Sintáxe: 

#
# Validação da entrada de dados na linha de comando
#
TAXON=$1	# Taxon path/filename or taxondir
BLASTDBNAME=$2	# Blastdb name
DBTYPE=$3	# Tipo de banco de dados nucl ou prot
if [[ $# -lt 3 ]]; then
	echo "Falta o caminho/nome ou o caminho do Taxon, o diretório do Blastdb a ser criado, ou o tipo do banco de dados!"
	echo "Sintáxe: ./fasta2blastdb.sh <TAXONDIR/TAXONFILENAME> <BLASTDBNAME> <BDTYPE: nucl/prot>"
	exit 0
fi

# Diretórios dos dados (Subjects)
# Salvar os arquivos contendo as sequencias referências no formato Fasta (obtidos do Genbank) 
# individualmente ou no formato multiseq Fasta neste diretório ou em sub-diretórios. Os arquivos
# serão contatenados recursivamente em um único arquivo refseq.fasta para criação do banco de dados
REFSEQDIR=${HOME}/data/REFSEQ   	# Habitualmente, podemos concentrar as REFSEQs aqui
# TAXDIR=${HOME}/data/REFSEQ/HEV   	# Se quisermos análisar o genoma do HEV apenas, por exemplo

# Diretório onde será criado o novo banco de dados refseq
BLASTDBDIR=${HOME}/data/BLASTDB/${BLASTDBNAME}
# BLASTDBDIR=${HOME}/data/BLASTDB/HEV_DB      	# # Se quisermos análisar o genoma do HEV apenas, por exemplo

# Reseta o diretório antes de criar um novo banco de dados
[[ -d ${BLASTDBDIR} ]] && rm -r ${BLASTDBDIR}
[[ ! -d ${BLASTDBDIR} ]] && mkdir -vp ${BLASTDBDIR}

# Remove o arquivo contendo as sequencias referência, se houver, antes de criar um novo
[[ -f ${REFSEQDIR}/${BLASTDBNAME}.fasta ]] && rm ${REFSEQDIR}/${BLASTDBNAME}.fasta

# Se TAXON for um diretório, concatena todos os arquivos .fasta em ${BLASTDBNAME}.fasta
# Exceto o arquivo ${BLASTDBNAME}.fasta que é gerado pelo make_refgen.sh
echo "Concatenando as sequencias referências em ${BLASTDBNAME}.fasta..."
if [ -f ${TAXON} ]; then
	cat ${TAXON} > "${BLASTDBDIR}/${BLASTDBNAME}.fasta"
else
	# find ${TAXON} -type f -iname '*.fasta' -print0 | sort -z | xargs -0 cat > "${BLASTDBDIR}/${BLASTDBNAME}.fasta"
	find ${TAXON} -name '*.fasta' -exec cat {} + > "${BLASTDBDIR}/${BLASTDBNAME}.fasta"
fi

# Processa a linha de descrição das sequencias referências para conter apenas o número de acesso sem espaços
echo "Processando os labels do arquivo ${BLASTDBNAME}.fasta..."
[[ -f ${BLASTDBDIR}/${BLASTDBNAME}.old ]] && rm ${BLASTDBDIR}/${BLASTDBNAME}.old
mv ${BLASTDBDIR}/${BLASTDBNAME}.fasta ${BLASTDBDIR}/${BLASTDBNAME}.old
while read -r line; do
	if echo "$line" | grep ">"; then
    		echo "$line" | cut -d "." -f 1 >> ${BLASTDBDIR}/${BLASTDBNAME}.fasta
	else
		echo "$line" >> ${BLASTDBDIR}/${BLASTDBNAME}.fasta
	fi
done < "${BLASTDBDIR}/${BLASTDBNAME}.old"

# Cria a lista de números de acc Genbank a partir do arquivo .fasta
echo "Criando o arquivo BLASTDBDIR.acc..."
[[ -f ${BLASTDBDIR}/${BLASTDBNAME}.acc ]] && rm  ${BLASTDBDIR}/${BLASTDBNAME}.acc
grep ">" ${BLASTDBDIR}/${BLASTDBNAME}.fasta | sed 's/>//' | cut -d " " -f 1 > ${BLASTDBDIR}/${BLASTDBNAME}.acc

# Cria a lista de taxid a partir nos números de acc Genbank
[[ -f ${BLASTDBDIR}/${BLASTDBNAME}.map ]] && rm  ${BLASTDBDIR}/${BLASTDBNAME}.map
# Retrive Taxid
 echo "Criando o arquivo ${BLASTDBNAME}.map..."
 while read -r line; do
 # echo "$line "$(efetch -db nuccore -id "$line" -format docsum | xtract -pattern DocumentSummary -element TaxId) >>${BLASTDBDIR}/${BLASTDBNAME}.map
	echo "$line "$(esearch -db assembly -q "$line" | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Taxid) >>${BLASTDBDIR}/${BLASTDBNAME}.map
done < ${BLASTDBDIR}/${BLASTDBNAME}.acc

# Alternativamente, podemos obter o Taxid usado esearch em combinação com esummary
# esearch -db assembly -q 'M62321.1' | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Taxid

# Cria o banco de dados refseq para busca pelos programas Blast a partir de um arquivo .fasta
echo "Criando o banco de dados data/blastdb/${BLASTDBNAME}..."
makeblastdb -in ${BLASTDBDIR}/${BLASTDBNAME}.fasta -parse_seqids -taxid_map ${BLASTDBDIR}/${BLASTDBNAME}.map -dbtype ${DBTYPE} -out ${BLASTDBDIR}/refseq
echo "Banco de dados criado com sucesso!"

# Faz o donwload do taxdb
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz -P ${BLASTDBDIR}

# IMPORTATE: Após incluir uma nova sequencia no diretório REFSEQ e executar uma nova análise, é importante
# atualizar os arquivos WIMP para sumarizar o relatório da busca pelo Blast. Abaixo seguem alguns reports
# para conferência
# echo "Conferir o arquivo refseq.fasta e os arquivos .wimp..."
# echo "Total de taxons encontrados no arquivo refseq.fasta atual: $(grep -c ">" ${REFSEQDIR}/refseq.fasta)"
# echo "Total de números de acesso (refseq.acc): $(wc -l ${REFSEQDIR}/refseq.acc)"
# echo "Total de números de acesso com taxid (refseq.map): $(wc -l ${REFSEQDIR}/refseq.map)"
# echo -e "Total de números de acesso em cada arquivo .wimp:\n $(wc -l ${HOME}/data/WIMP/*.wimp)"
# echo "IMPORTANTE: Atualizar os arquivos .wimp para uma relatoria correta!"
