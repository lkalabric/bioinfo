#!/bin/bash

# Cria o banco de dados Blast ou Diamond a partir de arquivos .fasta para busca
# Autor: Luciano Kalabric Silva
# Data: 24/03/2023
# Última atualização: 21/12/2023
# Log: Debugging

# Sintáxe: 

#
# Validação da entrada de dados na linha de comando
#
TAXON=$1	# Taxon path/filename or taxondir
DBNAME=$2	# Database name
DBTYPE=$3	# Tipo de banco de dados nucl ou prot
if [[ $# -lt 3 ]]; then
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

      # Reseta o diretório antes de criar um novo banco de dados
      [[ -d ${DBDIR} ]] && rm -r ${DBDIR}
      [[ ! -d ${DBDIR} ]] && mkdir -vp ${DBDIR}
      
      # Se TAXON for um diretório, concatena todos os arquivos .fasta em ${DBNAME}/refseq.fasta antes de montar o banco de dados
      echo "Concatenando as sequencias referências em ${DBDIR}/refseq.fasta..."
      if [ -f ${TAXON} ]; then
      	cat ${TAXON} > "${DBDIR}/refseq.fasta"
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
      [[ -f ${DBDIR}/refseq.map ]] && rm  ${DBDIR}/refseq.map
      # Retrive Taxid
       echo "Criando o arquivo ${DBDIR}/refseq.map..."
      while read -r line; do
      	echo "$line "$(esearch -db assembly -q "$line" < /dev/null | esummary | xtract -pattern DocumentSummary -element Taxid) >> ${DBDIR}/refseq.map
      done < ${DBDIR}/refseq.acc

      # Cria o banco de dados refseq propriamente dito para busca pelos programas Blast
      echo "Criando o banco de dados ${DBNAME}..."
      makeblastdb -in ${DBDIR}/refseq.fasta -parse_seqids -taxid_map ${DBDIR}/refseq.map -dbtype ${DBTYPE} -out ${DBDIR}/refseq
      echo "Banco de dados criado com sucesso!"

      exit 1
    ;;
    
    "-diamond")
      # Diretório onde será criado o novo banco de dados refseq
      DBDIR=${HOME}/data/DIAMONDDB/${DBNAME}
        
      # Reseta o diretório antes de criar um novo banco de dados
      [[ -d ${DBDIR} ]] && rm -r ${DBDIR}
      [[ ! -d ${DBDIR} ]] && mkdir -vp ${DBDIR}

      # Se TAXON for um diretório, concatena todos os arquivos .fasta em ${DBNAME}/refseq.fasta antes de montar o banco de dados
      echo "Concatenando as sequencias referências em ${DBDIR}/refseq.fasta..."
      if [ -f ${TAXON} ]; then
        cat ${TAXON} > "${DBDIR}/refseq.fasta"
      else
        # find ${TAXON} -type f -iname '*.fasta' -print0 | sort -z | xargs -0 cat > "${DBDIR}/refseq.fasta"
        find ${TAXON} -name '*.fasta' -exec cat {} + > "${DBDIR}/refseq.fasta"
      fi

      conda activate diamond
      # Cria o banco de dados refseq propriamente dito para busca pelos programas Blast
      echo "Criando o banco de dados ${DBNAME}..."
      diamond makedb --in ${DBDIR}/refseq.fasta --db refseq
      cp refseq.* ${DBDIR}/

      exit 2
    ;;
    
    *)
        echo "Invalid parameter!"
        exit 3
esac

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
