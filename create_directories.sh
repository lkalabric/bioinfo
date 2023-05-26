#!/bin/bash

# command name: create_directories.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 26 MAY 2023
# objetive: Create standardized directories to the Linux filesystem
# Syntax: ./create_directories.sh
# Link: Not available

# Suggestion of a directory tree
# bin/ - binary executable files
# data/ - your data files
# examples/ - example files
# repos/ - github repositories
# results/ - your results
# scripts/ - commands created in different scripting languages

# Since we need this command from bioinfo repository a priori, 
# we need to create a local repos/ directory using  mkdir  and 
# clone bioinfo using the following commands:
# $ mkdir repos/
# $ git clone https://github.com/lkalabric/bioinfo.git ./repos

# Command to create multiple directories at once in your HOME directory
mkdir ${HOME}/bin ${HOME}/data ${HOME}/examples ${HOME}/results ${HOME}/scripts
# Note: Feel free to include other directories of your choise
