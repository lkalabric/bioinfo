#!/bin/bash

#Testando o carregamento de um arquivo de configuração
CONFIG_FILE="./param/fastqc.param"
if [ -f "$CONFIG_FILE" ]; then
    # Carrega as variáveis do arquivo de configurações
    source "$CONFIG_FILE" # . ./param/fastqc.param 
else
    echo "ERRO: Arquivo de configuração '$CONFIG_FILE' não encontrado."
    exit 1
fi
echo "--- Carregamento de Variáveis Concluído ---"
echo "Filtro padrão: $FILTER"
echo "Número de threads personalizado: $THREADS"
