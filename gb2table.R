# 1. Instalação das dependências (descomente se necessário)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("genbankr")
  install.packages("tidyverse")

library(genbankr)
library(tidyverse)

# Função para extrair dados de um arquivo GenBank com múltiplos registros
get_gb_features <- function(gb_file_path) {
  
  # Lendo o arquivo GenBank (suporta arquivos com múltiplos LOCUS)
  gb_parsed <- readGenBank(gb_file_path, partial = TRUE)
  
  # Se o arquivo tiver múltiplos registros, readGenBank devolve uma lista
  if (!is(gb_parsed, "GenBankRecord")) {
    records <- gb_parsed
  } else {
    records <- list(gb_parsed)
  }
  
  # Iterar sobre cada registro e extrair os dados de CDS
  results <- map_df(records, function(rec) {
    acc <- locus(rec)@accession  # Accession Number
    cds_info <- cds(rec)         # Tabela de CDS (coding sequences)
    
    # Extrair colunas desejadas se existirem
    df <- as.data.frame(cds_info) %>%
      select(
        accession = any_of("accession"),
        protein_id = any_of("protein_id"),
        translation = any_of("translation")
      )
    
    # Preencher accession se não tiver vindo direto na tabela CDS
    if (!"accession" %in% colnames(df) || all(is.na(df$accession))) {
      df$accession <- acc
    }
    
    return(df)
  })
  
  return(results)
}

# --- Exemplo de uso ---
df_resultados <- get_gb_features("C:\Users\luciano.kalabric\Downloads\sequence.seq")
head(df_resultados)
# write.csv(df_resultados, "C:\Users\luciano.kalabric\Downloads\metadados.csv", row.names = FALSE)