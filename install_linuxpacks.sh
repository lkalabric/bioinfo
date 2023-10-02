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
# Function to test if package exists in Debian
function package_exists() {
    dpkg -s "$2" &> /dev/null
    return $?
}

# Linux packages are listed in a files *.packs at the following $PACKAGELIST_DIR
PACKAGELIST_DIR="${HOME}/repos/bioinfo"


# Validate parameters
if [ $# = 0 ]; then
	echo "Sintax: install_linuxpacks.sh <-i to install/-l to list> <package name or package list *.packs file>"
	exit 0;
else
	if [[ ! -f ${PACKAGELIST_DIR}/$2 ]]; then
		echo "Package name or package list *.packs file is required!"
		echo "Sintax: install_linuxpacks.sh <-i to install/-l to list> <package name or package list *.packs file>"
		exit 0
	else
 		PACKAGE_LIST = $(cat ${PACKAGELIST_DIR}/$2)
		case $1 in
			"-i" ) echo "Installation in progress..."
				# Pior to any installation it is recommended to update-upgrade your Linux Distro
				# Update & upgrade your Linux Distro
				echo "Updating & upgrading installed packages before starting any new installation..."
				sudo apt-get update
				sudo apt list --upgradable
				sudo apt-get upgrade
    			;;
			"-l" ) echo "Listing package(s) name(s) and descrition..."; for PACKAGE_NAME in "${PACKAGE_LIST[@]}"; do apt-cache search ^${PACKAGE_NAME}$; done; exit 0 ;;
			* ) echo "Invalid option!"; exit 0 ;;
		esac	
  		if [ -f ${PACKAGELIST_DIR}/$2 ]; then
			PACKAGE_LIST=($(cat ${PACKAGELIST_DIR}/$2))
		else
			PACKAGE_LIST=$2
			if ! package_exists ${PACKAGE_LIST}; then
				echo "Package name wrong or package list *.packs not found!"
				exit 0				
			else
				echo "Package ${PACKAGE_LIST} is available in the Debian Distro!"
			fi
		fi
	fi
fi

# Read package list and install each linux command if exists
for PACKAGE_NAME in "${PACKAGE_LIST[@]}"; do 	
	if ! which $PACKAGE_NAME > /dev/null; then
		echo -e "$PACKAGE_NAME is not installed! Install? (y/n) \c"
		read -r
		echo $REPLY
		if [[ $REPLY = "y" ]]; then
			sudo apt-get install ${PACKAGE_NAME}
			echo "`date` sudo apt-get install $PACKAGE_NAME" >> ${HOME}/logs/install_linuxpackages.log
			else
			echo "You can install it anytime!"
		fi
	else
		echo "$PACKAGE_NAME already installed in your Linux Distro!"
	fi
done
