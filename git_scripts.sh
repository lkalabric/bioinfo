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
		echo "Repository not present in repos/. Git clone it first!"
	else
		# Git pull REPO
  		echo "Git pulling ${REPO} repo..."
  		cd ${REPO_DIR}/${REPO}
		git pull
		# Reset scripts/ dir and copy files .sh and .R to it
  		rm -r ${SCRIPT_DIR}
    		mkdir ${SCRIPT_DIR}
  		find . \( -name '*.sh' -o -name '*.R' \) -exec cp {} ${SCRIPT_DIR} \;
		chmod +x ${SCRIPT_DIR}/*.sh
  		cd
    		export PATH="/${SCRIPT_DIR}:$PATH" Add to PATH Temporarily
      		# Better is to add Permanently, need to include the line above in the file .bashrc 
	fi	
 else
 	echo "List of cloned repositories:"
	ls $REPO_DIR
	# Git put all repos
	for dir in ${REPO_DIR}/*; do
 		echo "Git pulling ${dir} repo..."
  		cd ${dir}
		git pull
		# Reset scripts/ dir and copy files .sh and .R to it
  		rm -r ${SCRIPT_DIR}
    		mkdir ${SCRIPT_DIR}
  		find . \( -name '*.sh' -o -name '*.R' \) -exec cp {} ${SCRIPT_DIR} \;
		chmod +x ${SCRIPT_DIR}/*.sh
  		cd
    		export PATH="/${SCRIPT_DIR}:$PATH" # Add to PATH Temporarily
      		# Better is to add Permanently, need to include the line above in the file .bashrc 
	done
 fi
