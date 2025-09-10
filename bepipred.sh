#!/bin/bash

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

# Por convenção, 'exit 0' indica sucesso
exit 0
