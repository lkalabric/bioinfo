#!/bin/bash

# Cria o banco de dados Blast ou Diamond a partir de arquivos .fasta para busca
# Autor: Luciano Kalabric Silva
# Data: 24/03/2023
# Última atualização: 29/04/2024
# Log: Debugging
# 15/04/2024 - Permite que o script continue de onde parou para permitir concluir a criação do arquivo .map
# Sintáxe: Em revisão

#
# Validação da entrada de dados na linha de comando
#
DBTOOL=$1	# Feramneta de banco de dados Blast ou Diamond
TAXON=$2	# Taxon path/filename or taxondir
DBNAME=$3	# Database name
DBTYPE=$4	# Tipo de banco de dados nucl ou prot
if [[ $# -lt 4 ]]; then
	echo "Falta o caminho/nome ou o caminho do Taxon, o diretório do Blastdb a ser criado, ou o tipo do banco de dados!"
	echo "Sintáxe: ./fasta2db.sh <-blast/-diamond> <TAXONDIR/TAXONFILENAME> <DBNAME> <BDTYPE: nucl/prot>"
 	exit 0
fi

# Diretórios dos dados data/ (Subjects)
# Salvar os arquivos contendo as sequencias referências no formato Fasta (obtidos do Genbank) 
# individualmente ou no formato multiseq Fasta neste diretório ou em sub-diretórios. Os arquivos
# serão contatenados recursivamente em um único arquivo refseq.fasta para criação do banco de dados

case $1 in
    ##
    # Blast database
    ##
	"-blast")
		# Diretório onde será criado o novo banco de dados refseq
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
		
		# Se TAXON for um diretório, concatena todos os arquivos .fasta em ${DBDIR}/refseq.fasta antes de montar o banco de dados
		echo "Concatenando as sequencias referências ${TAXON} em ${DBDIR}/refseq.fasta..."
		if [ -f ${TAXON} ]; then
  			cp "${HOME}/${TAXON}" "${DBDIR}/refseq.fasta"
		else
			# find ${TAXON} -type f -iname '*.fasta' -print0 | sort -z | xargs -0 cat > "${DBDIR}/refseq.fasta"
			find ${TAXON} -name '*.fasta' -exec cat {} + > "${DBDIR}/refseq.fasta"
		fi
	
		# Processa a linha de descrição das sequencias referências para conter apenas o número de acesso sem espaços
		# echo "Processando os labels do arquivo ${DBDIR}/refseq.fasta..."
		# [[ -f ${DBDIR}/refseq.old ]] && rm ${DBDIR}/refseq.old
		# mv ${DBDIR}/refseq.fasta ${DBDIR}/refseq.old
		# while read -r line; do
		#	if echo "$line" | grep ">"; then
		#    		echo "$line" | cut -d "." -f 1 >> ${DBDIR}/refseq.fasta
		#	else
		#		echo "$line" >> ${DBDIR}/refseq.fasta
		#	fi
		# done < "${DBDIR}/refseq.old"

		# Cria a lista de números de acc Genbank a partir do arquivo .fasta
		echo "Criando o arquivo ${DBDIR}/refseq.acc..."
		[[ -f ${DBDIR}/refseq.acc ]] && rm  ${DBDIR}/refseq.acc
		grep ">" ${DBDIR}/refseq.fasta | sed 's/>//' | cut -d " " -f 1 > ${DBDIR}/refseq.acc
		
		# Cria a lista de taxid a partir nos números de acc Genbank
		echo "Criando o arquivo ${DBDIR}/refseq.map..."
  		[[ -f ${DBDIR}/refseq.map ]] && rm  ${DBDIR}/refseq.map
		touch ${DBDIR}/refseq.map
  		
    		# Retrive Taxid
		# Caso seja necessário continuar, armazenar o último acc processado
      		lastacc=$(tail -n 1 ${DBDIR}/refseq.map | cut -d " " -f 1)
		while read -r line; do
			# Caso seja necessário continuar, pula as linhas com os acc já processados
   			# [[ ! -z $(grep "$line" "${DBDIR}/refseq.map") ]] && continue
      			echo "Downloading taxid do acc $line..." | tee -a ${DBDIR}/refseq.log
			echo "$line "$(esearch -db assembly -query "$line" < /dev/null | esummary | xtract -pattern DocumentSummary -element Taxid) >> ${DBDIR}/refseq.map
		done < ${DBDIR}/refseq.acc
		
		# Cria o banco de dados refseq propriamente dito para busca pelos programas Blast
		echo "Criando o banco de dados ${DBNAME}..."
		makeblastdb -in ${DBDIR}/refseq.fasta -parse_seqids -taxid_map ${DBDIR}/refseq.map -dbtype ${DBTYPE} -out ${DBDIR}/refseq
		echo "Banco de dados criado com sucesso!"
		
		exit 1
	;;

"-diamond")
		# Validação do tipo de busca
    		if [ ${DBTYPE} == "-nucl" ]; then
			echo "Invalid parameter! DIAMOND is a sequence aligner for protein and translated DNA searches only."
   			exit 0
       		fi
		
  		# Diretório onde será criado o novo banco de dados reference
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
		
		# Se TAXON for um diretório, concatena todos os arquivos .fasta em ${DBDIR}/reference.fasta antes de montar o banco de dados
		echo "Concatenando as sequencias referências ${TAXON} em ${DBDIR}/refseq.fasta..."
		if [ -f ${TAXON} ]; then
  			cp "${HOME}/${TAXON}" "${DBDIR}/reference.fasta"
		else
			# find ${TAXON} -type f -iname '*.fasta' -print0 | sort -z | xargs -0 cat > "${DBDIR}/reference.fasta"
			find ${TAXON} -name '*.fasta' -exec cat {} + > "${DBDIR}/refseq.fasta"
		fi
	
		# Cria a lista de números de acc Genbank a partir do arquivo .fasta
		echo "Criando o arquivo ${DBDIR}/reference.acc..."
		[[ -f ${DBDIR}/reference.acc ]] && rm  ${DBDIR}/reference.acc
		grep ">" ${DBDIR}/reference.fasta | sed 's/>//' | cut -d " " -f 1 > ${DBDIR}/reference.acc
		
		# Cria a lista de taxid a partir nos números de acc Genbank
		echo "Criando o arquivo ${DBDIR}/reference.map..."
  		[[ -f ${DBDIR}/reference.map ]] && rm  ${DBDIR}/reference.map
		touch ${DBDIR}/reference.map
  		
    		# Retrive Taxid
		# Caso seja necessário continuar, armazenar o último acc processado
      		lastacc=$(tail -n 1 ${DBDIR}/reference.map | cut -d " " -f 1)
		while read -r line; do
			# Caso seja necessário continuar, pula as linhas com os acc já processados
   			# [[ ! -z $(grep "$line" "${DBDIR}/reference.map") ]] && continue
      			echo "Downloading taxid do acc $line..." | tee -a ${DBDIR}/reference.log
			echo "$line "$(esearch -db assembly -query "$line" < /dev/null | esummary | xtract -pattern DocumentSummary -element Taxid) >> ${DBDIR}/reference.map
		done < ${DBDIR}/reference.acc
		
		# Cria o banco de dados reference propriamente dito para busca pelos programas Blast
		echo "Criando o banco de dados ${DBNAME}..."
		diamond makedb -in ${DBDIR}/reference.fasta -d ${DBDIR}/reference
		echo "Banco de dados criado com sucesso!"
		
		exit 2
	;;		
      
*)
		echo "Invalid parameter!"
		exit 3
esac
