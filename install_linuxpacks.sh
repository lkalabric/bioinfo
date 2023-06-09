#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# URL:
# last update: 13 OUT 2021
# Objetive: Install Ubuntu packages and keep record of all installations
# Syntax: ./install_linuxpacks.sh <-i/-l> <package_name/package_list *.packs file>
# Link: https://stackoverflow.com/questions/1298066/how-can-i-check-if-a-package-is-installed-and-install-it-if-not

# This script is good to superuser or root user only!!!
if [[ $(sudo -v) ]]; then
    sudo -v
    exit 0
fi

# Linux packages listed in a file .packs
PACKAGELIST_DIR="${HOME}/repos/bioinfo"
if [[ -z $2 ]]; then
	echo "Package name or package list *.packs file is required!"
	echo "Syntax: ./install_linuxpacks.sh <-i/-l> <package name/package list *.packs file>"
	exit 0
else
	if [ "${2: -5}" == ".packs" ]; then
		if [ -f ${PACKAGELIST_DIR}/$2 ]; then
			PACKAGE_LIST=($(cat ${PACKAGELIST_DIR}/$2))
		else
		   echo "File $2 is not a .packs or does not exist."
		   exit 0
		fi
	else
		PACKAGE_LIST=$2
	fi
fi

echo $PACKAGE_LIST

# Validate parameters
if [ $# = 0 ]; then
	echo "Parameter wrong or missing!"
	echo "Sintax: install_linuxpacks.sh <-i/-l/-h/--help> <package name or list>"
	exit 0;
else
	case $1 in
		"-i" ) echo "Installation in progress..."; exit 0 ;;
		"-l" ) echo "Listing packages names and descrition..."; for PACKAGE_NAME in "${PACKAGE_LIST[@]}"; do apt-cache search ^${PACKAGE_NAME}$; done; exit 0 ;;
		* ) echo "Invalid option!"; exit 0 ;;
	esac
fi

# Read package list files and install each linux command if not exists
for PACKAGE_NAME in "${PACKAGE_LIST[@]}"; do 	
	if ! which $PACKAGE_NAME > /dev/null; then
		echo -e "$PACKAGE_NAME is not found! Install? (y/n) \c"
		read -r
		echo $REPLY
		if [[ $REPLY = "y" ]]; then
			sudo apt-get install ${PACKAGE_NAME}
			echo -ne "`date` sudo apt-get install $PACKAGE_NAME\r" >> ${HOME}/logs/install_linuxpackages.log
			else
			echo "You can install it anytime!"
		fi
		else
		echo -e "$PACKAGE_NAME already installed in your Linux Distro!"
	fi
done
