#!/bin/bash

# Author: Luciano Kalabric Silva
# Institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# Objective: Install and execute bepipred3
# Syntax: install_bepipred3.sh
# Link: https://github.com/UberClifford/BepiPred3.0-Predictor/tree/main

# Este script verifica se um ambiente Conda com um nome específico
# já existe. Se não existir, ele cria o ambiente.

# --- Configurações ---
# Defina o nome do ambiente Conda que você quer verificar/criar
# Exemplo: ENV_NAME="bepipred3"
ENV_NAME="bepipred3"

# Defina as bibliotecas que você quer instalar, se o ambiente for criado
# Exemplo: PYTHON_VERSION="python=3.9"
PYTHON_VERSION="python=3.8.8"
# Exemplo: PACKAGES="numpy pandas scikit-learn"
PACKAGES=""

# --- Lógica do Script ---
echo "Verificando a existência do ambiente Conda: '$ENV_NAME'..."

# O comando 'conda info --envs' lista todos os ambientes.
# O 'grep -q' busca o nome do ambiente e retorna 0 se encontrar, 1 se não.
# A opção '-q' (quiet) suprime a saída do grep.
if conda info --envs | grep -q "$ENV_NAME"; then
    echo "Sucesso: O ambiente '$ENV_NAME' já existe."
else
    echo "Ambiente Conda '$ENV_NAME' não encontrado. Criando..."

    # O comando 'conda create' cria o ambiente.
    # O '-y' (yes) aceita todas as confirmações automaticamente.
    if conda create --name "$ENV_NAME" $PYTHON_VERSION $PACKAGES -y; then
        echo "Sucesso: O ambiente '$ENV_NAME' foi criado com sucesso!"                     
    else
        echo "Erro: Falha ao criar o ambiente '$ENV_NAME'. Verifique as permissões ou a instalação do Conda."
        # Encerra o script com um código de erro
        exit 1
    fi
fi

REPO_URL="https://github.com/UberClifford/BepiPred3.0-Predictor.git"
REPO_DIR="~/repos/$ENV_NAME"
# Clonagem script para execução do bepipred3 do GitHub
echo "Clonando o repositório do bepipred3 localmente..."
if [ -d "$REPO_DIR" ]; then
    echo "Sucesso: O diretório '$REPO_DIR' já existe. Nenhuma ação de clone necessária."
else
    echo "Diretório '$DIR_NAME' não encontrado. Iniciando o processo de clonegem..."
    # O comando 'git clone' clona o repositório.
    if git clone "$REPO_URL"; then
        echo "Sucesso: Repositório clonado com sucesso para '$REPO_DIR'."
        echo "Atualizando o script bepipr3_CLI.py..."
        cp ~/repos/bepipred3/bepipred3_CLI.py ~/scripts/
        chmod +x ~/scripts/bepipred3_CLI.py
    else
        echo "Erro: Falha ao clonar o repositório. Verifique a URL e sua conexão com a internet."
        # Sair com um código de erro para indicar falha
        exit 1
    fi
fi

# Etapa final da instalação será manual, pois o comando conda activate não executa dentro do script!
echo "Para conclui a instalação, execute conda init e conda activate $ENV_NAME"
echo "Com o ambiente "$ENV_NAME" ativado, inicie a instalação do programa executável utilizando o comando:"
echo "pip3 install -r repos/bepipred3/requirements.txt"
        
# Por convenção, 'exit 0' indica sucesso
exit 0
