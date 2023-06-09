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
RD="${HOME}/repos/${REPO}"

# Scripts diretory
SD="${HOME}/scripts"

# Validate the parameter
if [ $# = 0 ]; then
	echo "Repository name required! Sintax: git_scripts.sh <repository>"
	echo "List of cloned repositories:"
	ls $RD
	exit 0;
else
	if [ ! -d ${RD} ]; then
		echo "Repository not present in repos/"
	else
		cd ${RD}
		git pull
		# Copy files only if they exist
  		find . \( -name '*.sh' -o -name '*.R' \) -exec cp {} ${SD} \;
		chmod +x ${SD}/*.sh
  		cd
	fi
fi
