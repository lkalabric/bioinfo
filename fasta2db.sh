#!/bin/bash

# Cria o banco de dados Blast ou Diamond a partir de arquivos .fasta para busca
# Autor: Luciano Kalabric Silva
# Data: 24/03/2023
# Última atualização: 29/04/2024
# Log: Debugging
# 15/04/2024 - Permite que o script continue de onde parou para permitir concluir a criação do arquivo .map
# Sintáxe: Em revisão
Pacotes Linux requeridos: ncbi-blast+ ou diamond-aligner
Pacotes de outros repositórios requeridos: Nenhum

#
# Validação da entrada de dados na linha de comando
#
DBTOOL=$1	# Ferramenta de banco de dados Blast ou Diamond
TAXON=$2	# Taxon path/filename or taxondir
DBNAME=$3	# Database name
DBTYPE=$4	# Tipo de banco de dados nucl ou prot
if [[ $# -lt 4 ]]; then
	echo "Falta o caminho/nome ou o caminho do Taxon, o diretório do Blastdb a ser criado, ou o tipo do banco de dados!"
	echo "Sintáxe: ./fasta2db.sh <DBTOOL: 'blast', 'diamond'> <TAXON: Path/File name.fasta> <DBNAME> <BDTYPE: 'nucl', 'prot'>"
 	exit 0
fi

# Diretórios dos dados data/ (Subjects)
# Salvar os arquivos contendo as sequencias referências no formato Fasta (obtidos do Genbank) 
# individualmente ou no formato multiseq Fasta neste diretório ou em subdiretórios. Os arquivos
# serão contatenados recursivamente em um único arquivo refseq.fasta para criação do banco de dados

case $1 in
    ##
    # Blast database
    ##
	"blast")
		# Diretório onde será criado o novo banco de dados blastdb
		DBDIR=${HOME}/data/BLASTDB/${DBNAME}
     		if [ -d ${DBDIR} ]; then
			read -n 1 -p "Diretório já existe, (R)esetar ou (C)ontinuar? " continuar
		else
   			mkdir -vp ${DBDIR}
		fi
  
		# Reseta o diretório antes de criar um novo banco de dados
		case $continuar in
		    	[Rr]) 
	      			echo -e "\nResetando o banco de dados..."
				rm -r ${DBDIR}
				mkdir -vp ${DBDIR}
      			;;
		    	[Cc]) 
       				echo -e "\nContinuando de onde paramos..."
       			;;
		esac
		
		# Se TAXON for um diretório, concatena todos os arquivos .fasta em ${DBDIR}/blastdb.fasta antes de montar o banco de dados
		echo "Concatenando as sequências referências ${TAXON} em ${DBDIR}/blastdb.fasta..."
		if [ -f ${TAXON} ]; then
  			cp "${HOME}/${TAXON}" "${DBDIR}/blastdb.fasta"
		else
			# find ${TAXON} -type f -iname '*.fasta' -print0 | sort -z | xargs -0 cat > "${DBDIR}/blastdb.fasta"
			find ${TAXON} -name '*.fasta' -exec cat {} + > "${DBDIR}/blastdb.fasta"
		fi
	
		# Processa a linha de descrição das sequencias referências para conter apenas o número de acesso sem espaços
		# echo "Processando os labels do arquivo ${DBDIR}/blastdb.fasta..."
		# [[ -f ${DBDIR}/blastdb.old ]] && rm ${DBDIR}/blastdb.old
		# mv ${DBDIR}/blastdb.fasta ${DBDIR}/blastdb.old
		# while read -r line; do
		#	if echo "$line" | grep ">"; then
		#    		echo "$line" | cut -d "." -f 1 >> ${DBDIR}/blastdb.fasta
		#	else
		#		echo "$line" >> ${DBDIR}/blastdb.fasta
		#	fi
		# done < "${DBDIR}/blastdb.old"

		# Cria a lista de números de acc Genbank a partir do arquivo .fasta
		echo "Criando o arquivo ${DBDIR}/blastdb.acc..."
		[[ -f ${DBDIR}/blastdb.acc ]] && rm  ${DBDIR}/blastdb.acc
		grep ">" ${DBDIR}/blastdb.fasta | sed 's/>//' | cut -d " " -f 1 > ${DBDIR}/blastdb.acc
		
		# Cria a lista de taxid a partir nos números de acc Genbank
		echo "Criando o arquivo ${DBDIR}/blastdb.map..."
  		[[ -f ${DBDIR}/blastdb.map ]] && rm  ${DBDIR}/blastdb.map
		touch ${DBDIR}/blastdb.map
  		
    		# Retrive Taxid
		# Caso seja necessário continuar, armazenar o último acc processado
      		lastacc=$(tail -n 1 ${DBDIR}/blastdb.map | cut -d " " -f 1)
		while read -r line; do
			# Caso seja necessário continuar, pula as linhas com os acc já processados
   			# [[ ! -z $(grep "$line" "${DBDIR}/blastdb.map") ]] && continue
      			echo "Downloading taxid do acc $line..." | tee -a ${DBDIR}/blastdb.log
			echo "$line "$(esearch -db assembly -query "$line" < /dev/null | esummary | xtract -pattern DocumentSummary -element Taxid) >> ${DBDIR}/blastdb.map
		done < ${DBDIR}/blastdb.acc
		
		# Cria o banco de dados blastdb propriamente dito para busca pelos programas Blast
		echo "Criando o banco de dados ${DBNAME}..."
		makeblastdb -in ${DBDIR}/blastdb.fasta -parse_seqids -taxid_map ${DBDIR}/blastdb.map -dbtype ${DBTYPE} -out ${DBDIR}/blastdb
		echo "Banco de dados criado com sucesso!"
		
		exit 1
	;;

"diamond")
		# Validação do tipo de busca
    		if [ ${DBTYPE} == "-nucl" ]; then
			echo "Invalid parameter! DIAMOND is a sequence aligner for protein and translated DNA searches only."
   			exit 0
       		fi
		
  		# Diretório onde será criado o novo banco de dados diamonddb
		DBDIR=${HOME}/data/DIAMONDDB/${DBNAME}

   		if [ -d ${DBDIR} ]; then
			read -n 1 -p "Diretório já existe, (R)esetar ou (C)ontinuar? " continuar
		else
   			mkdir -vp ${DBDIR}
       		fi
   		# Reseta o diretório antes de criar um novo banco de dados
		case $continuar in
		    	[Rr]) 
	      			echo -e "\nResetando o banco de dados..."
				rm -r ${DBDIR}
				mkdir -vp ${DBDIR}
      			;;
		    	[Cc]) 
       				echo -e "\nContinuando de onde paramos..."
       			;;
		esac
		
		# Se TAXON for um diretório, concatena todos os arquivos .fasta em ${DBDIR}/diamonddb.fasta antes de montar o banco de dados
		echo "Concatenando as sequências referências ${TAXON} em ${DBDIR}/refseq.fasta..."
		if [ -f ${TAXON} ]; then
  			cp "${HOME}/${TAXON}" "${DBDIR}/refseq.fasta"
		else
			# find ${TAXON} -type f -iname '*.fasta' -print0 | sort -z | xargs -0 cat > "${DBDIR}/refseq.fasta"
			find ${TAXON} -name '*.fasta' -exec cat {} + > "${DBDIR}/refseq.fasta"
		fi
	
		# Cria a lista de números de acc Genbank a partir do arquivo .fasta
		echo "Criando o arquivo ${DBDIR}/diamonddb.acc..."
		[[ -f ${DBDIR}/diamonddb.acc ]] && rm  ${DBDIR}/diamonddb.acc
		grep ">" ${DBDIR}/diamonddb.fasta | sed 's/>//' | cut -d " " -f 1 > ${DBDIR}/diamonddb.acc
		
		# Cria a lista de taxid a partir nos números de acc Genbank
		echo "Criando o arquivo ${DBDIR}/diamonddb.map..."
  		[[ -f ${DBDIR}/diamonddb.map ]] && rm  ${DBDIR}/diamonddb.map
		touch ${DBDIR}/diamonddb.map
  		
    		# Retrive Taxid
		# Caso seja necessário continuar, armazenar o último acc processado
      		lastacc=$(tail -n 1 ${DBDIR}/diamonddb.map | cut -d " " -f 1)
		while read -r line; do
			# Caso seja necessário continuar, pula as linhas com os acc já processados
   			# [[ ! -z $(grep "$line" "${DBDIR}/diamonddb.map") ]] && continue
      			echo "Downloading taxid do acc $line..." | tee -a ${DBDIR}/diamonddb.log
			echo "$line "$(esearch -db assembly -query "$line" < /dev/null | esummary | xtract -pattern DocumentSummary -element Taxid) >> ${DBDIR}/diamonddb.map
		done < ${DBDIR}/diamonddb.acc
		
		# Cria o banco de dados diamonddb propriamente dito para busca pelos programas Blast
		echo "Criando o banco de dados ${DBNAME}..."
		diamond makedb --in ${DBDIR}/diamonddb.fasta --db ${DBDIR}/diamonddb
		echo "Banco de dados criado com sucesso!"
		
		exit 2
	;;		
      
*)
		echo "Invalid parameter!"
		exit 3
esac
