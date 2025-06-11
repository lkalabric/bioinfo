#!/bin/bash

# Substitui os labels de arquivos multisequence FASTA
# Autor: Luciano Kalabric Silva (uso Gemini)
# Data: 11/06/2025
# Última atualização: 11/06/2025
# Log: Debugging
# Sintáxe: replace_FASTAlabels.sh arquivo_fasta_original.fasta novos_labels.txt

# Verifica se os argumentos necessários foram fornecidos
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Uso: $0 <arquivo_fasta_original> <arquivo_novos_labels>"
  echo "Exemplo: $0 sequencias.fasta novos_labels.txt"
  exit 1
fi

arquivo_fasta_original="$1"
arquivo_novos_labels="$2"
arquivo_fasta_modificado="$($1 > sed -e 's/\.fasta$//')_modificado.fasta" # Nome de arquivo temporário com timestamp

# Verifica se os arquivos existem e são legí­veis
if [ ! -f "$arquivo_fasta_original" ] || [ ! -r "$arquivo_fasta_original" ]; then
  echo "Erro: O arquivo FASTA original '$arquivo_fasta_original' não existe ou nÃ£o pode ser lido."
  exit 1
fi

if [ ! -f "$arquivo_novos_labels" ] || [ ! -r "$arquivo_novos_labels" ]; then
  echo "Erro: O arquivo de novos labels '$arquivo_novos_labels' não existe ou nÃ£o pode ser lido."
  exit 1
fi

echo "Iniciando a substituição de labels em '$arquivo_fasta_original' usando '$arquivo_novos_labels'..."

# Ler os novos labels em um array
# readarray é específico do Bash 4+; se você tiver Bash 3, use IFS=$'\n' read -d '' -r -a novos_labels < "$arquivo_novos_labels"
mapfile -t novos_labels < "$arquivo_novos_labels"

# Contador para o índice dos novos labels
label_index=0

# Processar o arquivo FASTA
# Usa 'while read line' para ler o arquivo linha por linha
while IFS= read -r line; do
  # Se a linha começa com '>', há um cabeçalho de sequência
  if [[ "$line" == ">"* ]]; then
    # Verifica se ainda há novos labels disponí­veis
    if [ "$label_index" -lt "${#novos_labels[@]}" ]; then
      # Imprime o novo label prefixado com '>'
      echo ">${novos_labels[$label_index]}" >> "$arquivo_fasta_modificado"
      label_index=$((label_index + 1))
    else
      # Se não houver mais novos labels, usa o label original e imprime um aviso
      echo "$line" >> "$arquivo_fasta_modificado"
      echo "Aviso: Não há mais novos labels disponíveis para o label original: $line. Mantendo o original." >&2
    fi
  else
    # Se não for um cabeçalho, imprime a linha da sequência como está
    echo "$line" >> "$arquivo_fasta_modificado"
  fi
done < "$arquivo_fasta_original"

echo "---------------------------------------------------------"
echo "Substituição concluí­da. O arquivo modificado foi salvo como:"
echo "$arquivo_fasta_modificado"
echo "---------------------------------------------------------"

# Você pode renomear o arquivo modificado para o nome original, se desejar
# mv "$arquivo_fasta_modificado" "$arquivo_fasta_original"
# echo "Arquivo original substituÃ­do por '$arquivo_fasta_modificado'."
