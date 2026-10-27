#!/bin/bash

# autor: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# URL:
# last update: 27 OUT 2026
# Objective: Updade/Upgrade and install Linux packages, and keep a record of all installations
# Syntax: ./install_linuxpacks.sh <-u/-i/-l> <package_name/package_list *.packs file>
# Link: https://stackoverflow.com/questions/1298066/how-can-i-check-if-a-package-is-installed-and-install-it-if-not

# Linux packages are listed in files *.packs in the packages files dir
PACKAGE_DIR="${HOME}/repos/bioinfo/"

# Cria o diretório de log, se ele não existe
mkdir -p "${HOME}/logs

# This script is for superuser or root user only!!!
#if [ "$EUID" -ne 0 ]
#  then echo "Please run as root"
#  exit
#fi
if ((EUID != 1000)); then
	echo "Please run as root"
	exit 0
fi
sudo -v

# Function to test if a package is installed on your Debian machine
function is_installed() {
     dpkg --verify "$1" 2>/dev/null
     return $?
}

# Validate parameters
if [ $# = 0 ]; then
	echo "Syntax: install_linuxpacks.sh <-u to update/upgrade / -i to install packages from command line / -l to install packages from a package list *.packs file>"
	exit 0;
else
	case $1 in
		"-u" ) echo "Update/Upgrade in progress..."
			# Prior to any installation, it is recommended to update/upgrade your Linux distro.
			echo "Updating & upgrading installed packages before starting any new installation..."
			sudo apt-get update
			sudo apt list --upgradable
			sudo apt-get upgrade 
			;;
		"-i" ) PACKAGE_LIST=$2
			echo "Installation in progress from package list ${PACKAGE_LIST}..."
			for PACKAGE_NAME in ${PACKAGE_LIST}; do
				# Check if package is installed and install it if not
				if ! is_installed "${PACKAGE_NAME}"; then
					echo -e "${PACKAGE_NAME} is not installed! Install? (y/n) \c"
					read -r
					echo $REPLY
					if [[ $REPLY == [Yy] ]]; then
						sudo apt-get install "${PACKAGE_NAME}"
						echo "`date` sudo apt-get install $PACKAGE_NAME" >> ${HOME}/logs/install_linuxpackages.log
					else
						echo "You can install it anytime!"
					fi
				else
					echo "${PACKAGE_NAME} already installed in your Linux Distro!"
				fi
			done
			;;
		"-l" ) PACKAGELIST_FILENAME=$2
			echo "Installation in progress from package file ${PACKAGELIST_FILENAME}..."
			while IFS= read -r PACKAGE_NAME || [[ -n "${$PACKAGE_NAME}" ]]; do
				# Pula linhas vazias ou que começam com #
				[[ -z "${$PACKAGE_NAME}" || "${$PACKAGE_NAME}" =~ ^# ]] && continue
				
				# Remove espaços em branco extras, se houver
				PACKAGE_NAME=$(echo "${$PACKAGE_NAME}" | xargs)

				# O restante da sua lógica de verificação e instalação aqui
				if ! is_installed "${$PACKAGE_NAME}"; then
					echo -e "${$PACKAGE_NAME} is not installed! Install? (y/n) \c"
					read -r
					echo $REPLY
					if [[ $REPLY == [Yy] ]]; then
						sudo apt-get install "${$PACKAGE_NAME}"
						mkdir -p logs # Cria este diretório se ele não existe
						echo "`date` sudo apt-get install $PACKAGE_NAME" >> ${HOME}/logs/install_linuxpackages.log
					else
						echo "You can install it anytime!"
					fi
				else
					echo "$PACKAGE_NAME already installed in your Linux Distro!"
				fi
			done < "${PACKAGE_DIR}/${PACKAGELIST_FILENAME}"
			;;
		* ) echo "Invalid option!"; exit 0;;
	esac	
fi
