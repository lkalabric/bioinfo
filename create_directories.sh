#!/bin/bash

# command name: create_directories.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 26 MAY 2023
# objetive: Create standardized directories and prepare your Linux filesystem
# Syntax: ./create_directories.sh
# Link: https://linuxize.com/post/how-to-add-directory-to-path-in-linux/

# Suggestion of a directory tree
# bin/ - binary executable files
# data/ - your data files
# examples/ - example files
# logs/ - log files
# repos/ - github repositories
# results/ - your results
# scripts/ - commands created in different scripting languages

# Since we need this command from bioinfo repository a priori, 
# we need to create a local repos/ directory using  mkdir  and 
# clone bioinfo using the following commands:
# $ mkdir repos/
# $ git clone https://github.com/lkalabric/bioinfo.git ./repos

# Command to create multiple directories at once in your HOME directory
mkdir ${HOME}/bin ${HOME}/data ${HOME}/examples ${HOME}/logs ${HOME}/results ${HOME}/scripts
# Note: Feel free to include other directories of your choise

# Adding a directory to PATH 
# The most common directories that hold executable programs are /bin, /sbin, /usr/sbin, /usr/local/bin and /usr/local/sbin
# To add scripts/ in PATH use the following command
# $ export PATH="$HOME/bin:$PATH"
# To make the change permanent, you need to define the $PATH variable in the shell configuration files
# $ nano ~/.bashrc # and append the following lines:
# # add scripts/ to PATH
# export PATH="$HOME/bin:$PATH"
# After save the changes, run:
# $ source ~/.bashrc

# Once all directories have been created one can install
# our scritps into the scripts/ directory running the following:
# cd
# $ bash repos/bioinfo/git_scripts.sh
