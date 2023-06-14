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
RD="${HOME}/repos"

# Scripts diretory
SD="${HOME}/scripts"

# Validate the parameter
if [ $# = 0 ]; then
	echo "Repository name required! Sintax: git_scripts.sh <repository>"
	echo "List of cloned repositories:"
	ls $RD
	exit 0;
else
	if [ ! -d "${RD}/${REPO}" ]; then
		echo "Repository not present in repos/"
	else
		cd ${RD}/${REPO}
		git pull
		# Copy files only if they exist
  		find ${RD}/${REPO} -name '*.sh' -exec cp {} ${SD} \;
		chmod +x ${SD}/*.sh
  		# Copy files only if they exist
    		find ${RD}/${REPO} -name '*.R' -exec cp {} ${SD} \;
      		cd
	fi
fi
