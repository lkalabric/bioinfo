#!/bin/bash

# autor: Luciano Kalabric Silva
# creation date (YYYY-MM-DD): 2023-05-08
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# links:


# Objetive: Organize files and diretories
# Syntax: ./organize_filesystem.sh

# Last update (YYYY-MM-DD): 13 OUT 2021
# Version control

# Declare an array with the directory's names to be created
declare -a DIRECTORY_LIST=('bin' 'data' 'examples' 'logs' 'repos' 'results' scripts')

# Change to home directory
cd ${HOME}

# Read each directory name from the array and try to create a directory
for DIRECTORY_NAME in ${DIRECTORY_LIST}; do 
  [ -d $DIRECTORY_NAME ] && echo "Directory ${$DIRECTORY_NAME} exists!" || mkdir $DIRECTORY_NAME
done


