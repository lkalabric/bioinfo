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
WD="${HOME}/repos"

# Scripts diretory
SD="${HOME}/scripts"

# Validate the parameter
if [ $# = 0 ]; then
    echo "Repository name required! Sintax: git_scripts.sh <repository>"
    echo "List of cloned repositories:"
    ls $WD
    exit 0;
else
	if [ ! -d "${WD}/${REPO}" ]; then
	    echo "Repository not present in repos/"
	else
	    cd ${WD}/${REPO}
	    git pull
	    cp *.sh *.R ${SD}
	    chmod +x ${SD}/*.sh
	    cd
     fi
fi
