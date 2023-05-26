#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 26 MAY 2023
# Objetive: Update scripts from a git repo
# Syntax: ./git_scripts.sh <repo>
# Link: 

# Passes repo name to the script
REPO=$1

# Validate the parameter
if [ $# = 0 ]; then
    echo "Repository name required! Sintax: git_scripts.sh <repository>"  
    exit 0;
else
	if [ ! -d "${HOME}/repos/${REPO}" ]; then
	    echo "Repository not present in repos/"
	else
	    cd ${HOME}/repos/${REPO}
	    git pull
	    cp *.sh ~/scripts/
	    chmod +x ~/scripts/*
	    cd
	fi
fi
