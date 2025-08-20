#!/bin/bash

# Este script realiza o mapeamento de sequências de peptídeos em sequencias de proteínas utilizando o app peptidemapper
EXAMPLES_DIR="examples/peptidemapper/"


wget http://genesis.ugent.be/maven2/com/compomics/utilities/5.0.39/utilities-5.0.39.zip
unzip utilities-5.0.39.zip
cd utilities-5.0.39
#java -cp utilities-5.0.39.jar com.compomics.cli.peptide_mapper.PeptideMapperCLI -p exampleFiles/PeptideMapping/yeast.fasta exampleFiles/PeptideMapping/yeast-pep-1k.csv results.csv
java -cp utilities-5.0.39.jar com.compomics.cli.peptide_mapper.PeptideMapperCLI -p $EXAMPLES_DIR/PeptideMapping/yeast.fasta $EXAMPLES_DIR/PeptideMapping/yeast-pep-1k.csv results.csv
