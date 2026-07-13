#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 26 MAY 2023
# Objetive: Update scripts from a git repo
# Syntax: ./git_scripts.sh <repo>
# Link: 

# Passes repo name to the script
REPO=$1

# Repository directory
REPO_DIR="${HOME}/repos"

# Scripts directory
SCRIPT_DIR="${HOME}/scripts"

# Validate the parameter
if [ $# = 1 ]; then
	if [ ! -d ${REPO_DIR}/${REPO} ]; then
		echo "Repository not present in repos/. Git clone it first!"
	else
		# Git pull REPO
  		echo "Git pulling ${REPO} repo..."
  		cd ${REPO_DIR}/${REPO}
		git pull		
	fi	
else
 	echo "List of cloned repositories:"
	ls ${REPO_DIR}
 	cd ${REPO_DIR}
 	# Git put all repos
	for dir in $(find . -mindepth 1 -maxdepth 1 -type d); do
 		REPO="${dir#./}"     
 		echo "Git pulling ${REPO} repo..."
  		cd ${REPO}
		git pull		
  		cd ..
	done
fi

# Garante que a pasta física realmente exista
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Criando o diretório $SCRIPT_DIR..."
    mkdir -p "$SCRIPT_DIR"
fi

# Define o arquivo de configuração do Bash
ARQUIVO_CONFIG="$HOME/.bashrc"

# Comando de export que queremos adicionar
LINHA_EXPORT="export PATH=\"\$PATH:$SCRIPT_DIR\""

# Verifica se o caminho já está presente no arquivo para evitar duplicatas
if grep -Fq "$SCRIPT_DIR" "$ARQUIVO_CONFIG"; then
    echo "O caminho '$SCRIPT_DIR' já está configurado no seu $ARQUIVO_CONFIG."
else
    # Adiciona a linha ao final do arquivo de configuração
    echo -e "\n# Adicionado pelo script de configuracao de PATH\n$LINHA_EXPORT" >> "$ARQUIVO_CONFIG"
    echo "Sucesso! O caminho 'SCRIPT_DIR' foi adicionado ao seu $ARQUIVO_CONFIG."
    echo "Para aplicar as alterações imediatamente, execute: source $ARQUIVO_CONFIG"
fi

# Reset scripts/ dir and copy files .sh and .R to it
echo "Restoring scripts and privileges..."
find "${REPO_DIR}/" -maxdepth 2 \( -name '*.sh' -o -name '*.R' -o -name '*.py' \) -exec cp {} ${SCRIPT_DIR} \;
find "${SCRIPT_DIR}/" -maxdepth 2 \( -name '*.sh' -o -name '*.R' -o -name '*.py' \) -exec chmod +x {} \;
