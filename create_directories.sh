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

# Since we need to clone at least bioinfo repo from github
# We need to create and clone this command a priori using
# $ mkdir repos
# $ git clone 

# Command to create multiple directories at once in your HOME directory
mkdir ${HOME}/bin ${HOME}/data ${HOME}/examples ${HOME}/repos ${HOME}/results ${HOME}/scripts
