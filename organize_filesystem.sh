#!/bin/bash

# Autor: Luciano Kalabric Silva
# Creation date (YYYY-MM-DD): 2023-05-08
# Institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# Objetive: Organize files and diretories
# Syntax: ./organize_filesystem.sh

# Last update (YYYY-MM-DD): 13 OUT 2021
# Version control

# Declare an array with the directory's names to be created
declare -a DIRECTORY_LIST=('bin' 'data' 'examples' 'logs' 'repos' 'results' 'scripts')

# Change to home directory
cd ${HOME}

# Read each directory name from the array and try to create a directory
for DIRECTORY_NAME in ${DIRECTORY_LIST[@]}; do 
  [ -d ${DIRECTORY_NAME} ] && echo "Directory ${$DIRECTORY_NAME} exists!" || mkdir ${DIRECTORY_NAME}
done

# Adding a directory to your PATH (temporary solution)
# Link: https://linuxize.com/post/how-to-add-directory-to-path-in-linux/
# export PATH="$HOME/bin:$PATH"
# export PATH="$HOME/scripts:$PATH"

# Add directories to PATH (permanent change)
if [-z $(grep -q organize_filesystem.sh .bashrc)]; then 
   echo "PATH already configured!"
else
  echo -e "\n" >> ${HOME}/.bashrc
  echo "# New directories in PATH included by script organize_filesystem.sh" >> ${HOME}/.bashrc
  echo "PATH="$HOME/bin:$HOME/scripts:$PATH"" >> ${HOME}/.bashrc
  echo "PATH configured successfully!"
  # Reload .bashrc without need to wait until next login
  source ${HOME}/.bashrc
fi


# IMPORTANT: Only run those line once, otherwise they will be repeated in the .bashrc file!

