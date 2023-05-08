#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# URL:
# last update: 13 OUT 2021
# Objetive: Install Ubuntu packages and keep record of all installations
# Syntax: ./install_linuxpackages.sh <command>
# Link: https://stackoverflow.com/questions/1298066/how-can-i-check-if-a-package-is-installed-and-install-it-if-not

COMMAND=$1
HOME=/home/kalabric

if [ $# = 0 ]; then
	echo "Sintax error: ./install_linuxpackages.sh <command>"
	exit 0
fi

if ! which $COMMAND > /dev/null; then
	echo -e "$COMMAND is not found! Install? (y/n) \c"
	read $REPLY
	if [ $REPLY = "y" ]; then
		sudo apt-get install ${COMMAND}
		echo -ne "`date` sudo apt-get install $COMMAND\r" >> ${HOME}/scripts/install_linuxpackages.log
	fi
fi
