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
REPO_DIR=${HOME}/repos

# Scripts directory
SCRIPT_DIR="${HOME}/scripts"

# Validate the parameter
if [ $# = 1 ]; then
	if [ ! -d ${REPO_DIR}/${REPO} ]; then
		echo "Repository not present in repos/"
	else
		# Git pull REPO
  		cd ${REPO_DIR}/${REPO}
		git pull
		# Copy files only if they exist
  		find . \( -name '*.sh' -o -name '*.R' \) -exec cp {} ${SCRIPT_DIR} \;
		chmod +x ${SCRIPT_DIR}/*.sh
  		cd
	fi	
 else
 	echo "List of cloned repositories:"
	ls $REPO_DIR
	# Git put all repos
	for dir in ${REPO_DIR}/*; do
 		echo "Git pulling ${dir}..."
  		cd ${dir}
		git pull
		# Copy files only if they exist
  		find . \( -name '*.sh' -o -name '*.R' \) -exec cp {} ${SCRIPT_DIR} \;
		chmod +x ${SCRIPT_DIR}/*.sh
  		cd
	done
 fi
