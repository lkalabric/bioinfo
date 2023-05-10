#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# URL:
# last update: 13 OUT 2021
# Objetive: Install Ubuntu packages and keep record of all installations
# Syntax: ./install_linuxpackages.sh <package>
# Link: https://stackoverflow.com/questions/1298066/how-can-i-check-if-a-package-is-installed-and-install-it-if-not

# Linux packages listed in a file .packs

# Validate parameters
case $1 in
	"--help") echo "Sintax: ./install_linuxpackages.sh <-i/-a/-h/--help> <package_list .packs>" ;;
	"-h") echo "Sintax: ./install_linuxpackages.sh <-i/-a/-h/--help> <package_list .packs>" ;;
	"-i") echo "Installation in progress..." ;;
	"-a") echo "Listing packages names and descrition..." ;;
	*) echo "Invalid option!" ;;
esac

# If no parameter, assume all! 
if [ $# = 0 ]; then
	$PACKAGE_LIST="${HOME}/repos/biolinux/*.packs"
	else
	PACKAGE_LIST=$1
fi

# Read package list(s) and install linux command if not exists
while IFS= read -r PACKAGE; do 
	echo $PACKAGE
	if ! which $PACKAGE > /dev/null; then
		echo -e "$PACKAGE is not found! Installation in progress..."
		# echo -e "$PACKAGE is not found! Install? (y/n) \c"
		# read -r
		# echo $REPLY
		#if [[ $REPLY = "y" ]]; then
			sudo apt-get install ${PACKAGE}
			echo -ne "`date` sudo apt-get install $PACKAGE\r" >> ${HOME}/logs/install_${PACKAGE}.log
		#	else
		#	echo "You can install it anytime!"
		#fi
		else
		echo -e "$PACKAGE already installed in your Linux Distro!"
	fi
done < ${PACKAGE_LIST}
