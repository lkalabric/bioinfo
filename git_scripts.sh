#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 26 MAY 2023
# Objetive: Update scripts from a git repo
# Syntax: ./git_scripts.sh <repo>
# Link: 

# Passes repo name to the script
REPO=$1

# Repository directory
REPO_DIR="${HOME}/repos"

# Scripts directory
SCRIPT_DIR="${HOME}/scripts"

# Validate the parameter
if [ $# = 1 ]; then
	if [ ! -d ${REPO_DIR}/${REPO} ]; then
		echo "Repository not present in repos/. Git clone it first!"
	else
		# Git pull REPO
  		echo "Git pulling ${REPO} repo..."
  		cd ${REPO_DIR}/${REPO}
		git pull		
	fi	
else
 	echo "List of cloned repositories:"
	ls ${REPO_DIR}
 	cd ${REPO_DIR}
 	# Git put all repos
	for dir in $(find . -mindepth 1 -maxdepth 1 -type d); do
 		REPO="${dir#./}"     
 		echo "Git pulling ${REPO} repo..."
  		cd ${REPO}
		git pull		
  		cd ..
	done
fi
# Reset scripts/ dir and copy files .sh and .R to it
echo "Restoring scripts and privileges..."
rm -r ${SCRIPT_DIR}
mkdir ${SCRIPT_DIR}
find "${REPO_DIR}/" -maxdepth 2 \( -name '*.sh' -o -name '*.R' -o -name '*.py' \) -exec cp {} ${SCRIPT_DIR} \;
find "${SCRIPT_DIR}/" -maxdepth 2 \( -name '*.sh' -o -name '*.R' -o -name '*.py' \) -exec chmod +x {} \;
# To add SCRIPT_DIR permanently in the PATH, one needs to edit .bashrc file and add the following line
# export PATH="/${SCRIPT_DIR}:$PATH"

