#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# URL:
# last update: 13 OUT 2021
# Objetive: Install Ubuntu packages and keep record of all installations
# Syntax: ./install_linuxpackages.sh <package>
# Link: https://stackoverflow.com/questions/1298066/how-can-i-check-if-a-package-is-installed-and-install-it-if-not

# Linux packages listed in a file .packs
PACKAGE_LIST=$1

# Validate parameters
if [ $# = 0 ]; then
	echo "Sintax error: ./install_linuxpackages.sh <package_list .packs>"
	exit 0
fi

# Read package list and install linux command if not exists
while IFS= read -r PACKAGE; do 
	if ! which $PACKAGE > /dev/null; then
		echo -e "${PACKAGE} is not found! Install? (y/n) \c\n"
		read REPLY
		if [ $REPLY = "y" ]; then
			sudo apt-get install ${PACKAGE}
			echo -ne "`date` sudo apt-get install $PACKAGE\r" >> ${HOME}/logs/install_${PACKAGE}.log
		fi
		else
		echo -e "$PACKAGE already installed in your Linux Distro!"
	fi
done < ${PACKAGE_LIST}
