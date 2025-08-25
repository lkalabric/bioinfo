#!/bin/bash

# command name: create_directories.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 26 MAY 2023
# objetive: Create standardized directories and prepare your Linux filesystem
# Syntax: ./create_directories.sh
# Link: https://linuxize.com/post/how-to-add-directory-to-path-in-linux/

# Suggestion of a directory tree
# bin/ - binary executable files
# data/ - your data files
# examples/ - example files
# logs/ - log files
# repos/ - github repositories
# results/ - your results
# scripts/ - commands created in different scripting languages
# scratch/ - temporary files
DIRECTORY_LIST=("bin" "data" "examples" "logs" "repos" "results" "scripts" "scratch")

# Check if the directory list is valid
if [ ! -f "$DIRECTORY_LIST" ] || [ ! -r "$DIRECTORY_LIST" ]; then
    echo "Erro: A lista '$DIRECTORY_LIST' não existe ou não é válida!"
    exit 1
fi
# Lê a variável contendo a lista de diretórios e cria cada diretório
echo "Iniciando a criação dos diretórios..."
while IFS= read -r DIR_NOME; do
    # O comando 'mkdir -p' cria o diretório (e todos os pais, se necessário).
    # A opção '-p' também evita erros se o diretório já existir.
    if mkdir -p "$DIR_NOME"; then
        echo "Sucesso: Diretório '$DIR_NOME' criado."
    else
        echo "Erro: Falha ao criar o diretório '$DIR_NOME'."
    fi
done < "$DIRECTORY_LIST"
echo "Processo concluído."
exit 0

# Since we need this command from bioinfo repository a priori, 
# we need to clone it into a local repos/bioinfo using the following commands:
# $ git clone https://github.com/lkalabric/bioinfo.git repos/bioinfo

# Now, let´s do the magic:
# $ bash repos/bioinfo/create_directories.sh

# Command to create multiple directories at once in your HOME directory
mkdir ${HOME}/bin ${HOME}/data ${HOME}/examples ${HOME}/logs ${HOME}/results ${HOME}/scripts
# Note: Feel free to include other directories of your choise

# Adding a directory to PATH 
# The most common directories that hold executable programs are /bin, /sbin, /usr/sbin, /usr/local/bin and /usr/local/sbin
# To add scripts/ in PATH use the following command
# $ export PATH="$HOME/scripts:$PATH"
# To make the change permanent, you need to define the $PATH variable in the shell configuration files
# $ nano ~/.bashrc # and append the following lines:
# # add scripts/ to PATH
# export PATH="$HOME/scripts:$PATH"
# After save the changes, run:
# $ source ~/.bashrc

# Once all directories have been created one can install
# our scritps into the scripts/ directory running the following:
# cd
# $ bash repos/bioinfo/git_scripts.sh
