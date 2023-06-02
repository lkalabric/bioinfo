#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# URL:
# last update: 13 OUT 2021
# Objetive: Install Ubuntu packages and keep record of all installations
# Syntax: ./install_linuxpacks.sh <-i/-a/-h/--help> <package_name/package_list.packs/*>
# Link: https://stackoverflow.com/questions/1298066/how-can-i-check-if-a-package-is-installed-and-install-it-if-not

# This script is good to superuser or root user only!!!
# Need to include an error message in case the user does not have privilegies!!!
if [[ $(sudo -v) ]]; then
    sudo -v
    exit 0
fi


# Linux packages listed in a file .packs
PACKAGELIST_DIR="${HOME}/repos/bioinfo"

# Validate parameters
if [ $# = 0 ]; then
	echo "Package name (single installation)  or package list filename (batch installation) required! Sintax: install_linuxpacks.sh <package name or list>"  
	exit 0;
else
	if [ $# = 1 ]; then
		case $1 in
			"--help") echo "Sintax: ./install_linuxpackages.sh <-i/-a/-h/--help> <filename.packs>"; exit 0 ;;
			"-h") echo "Sintax: ./install_linuxpackages.sh <-i/-a/-h/--help> <filename.packs>"; exit 0 ;;
			"-i") echo "Installation in progress..."; exit 0 ;;
			"-a") echo "Listing packages names and descrition..."; exit 0 ;;
		#	*) echo "Invalid option!"; exit 0 ;;
		esac
	else
		if [ "$2" == "*.packs" ]; then
        		PACKAGE_LIST=($(cat ${HOME}/repos/bioinfo/$2))
		else
			PACKAGE_LIST=$2
		fi
	fi
fi

# Read package list files and install each linux command if not exists
for PACKAGE_NAME in "${PACKAGE_LIST[@]}"; do 	
	if ! which $PACKAGE_NAME > /dev/null; then
		echo -e "$PACKAGE_NAME is not found! Installation in progress..."
		# echo -e "$PACKAGE_NAME is not found! Install? (y/n) \c"
		# read -r
		# echo $REPLY
		#if [[ $REPLY = "y" ]]; then
			sudo apt-get install ${PACKAGE_NAME}
			echo -ne "`date` sudo apt-get install $PACKAGE_NAME\r" >> ${HOME}/logs/install_linuxpackages.log
		#	else
		#	echo "You can install it anytime!"
		#fi
		else
		echo -e "$PACKAGE_NAME already installed in your Linux Distro!"
	fi
done
