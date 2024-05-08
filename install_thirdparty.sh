#!/bin/bash

# command name: install_thirdparty.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 13 JUN 2023
# objetive: Install thirdparty commands in Linux Ubuntu
# Syntax: ./install_thirdparty.sh

echo "Tem erro aqui?!"

# Update & upgrade your Linux Distro
# This script is good to superuser or root user only!!!
if [[ $(sudo -v) ]]; then
	echo "Updating & upgrading installed packages before starting any new installation..."
	sudo apt-get update
	sudo apt list --upgradable
	sudo apt-get upgrade
fi

# Install Miniconda, if not present
if [ ! -f ~/Downloads/Miniconda3-latest-Linux-x86_64.sh ]; then
	# Install and configure Miniconda if not installed
 	# Link: https://www.cyberithub.com/how-to-install-miniconda-on-ubuntu-20-04-lts-focal-fossa/
	echo "Downloading and installing Conda in your system..."
 	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -P Downloads/
	chmod +x Downloads/Miniconda3-latest-Linux-x86_64.sh
	bash Downloads/Miniconda3-latest-Linux-x86_64.sh
	# Configure PATH
	# To make the change permanent, you need to define the $PATH variable in a shell configuration file like ~/.bashrc.
	echo 'export PATH="$HOME/miniconda3/bin:$PATH" # add miniconda3/bin to PATH' >> ~/.bashrc # appends the export to ~/.bashrc file
	# After saving the file, run the following command to the export take effect:
	source ~/.bashrc
fi

# Using Miniconda
# To keep compatibility between Python versions, we strongly recommend creating different envs for each app
# Installation Bioconda packages
# https://anaconda.org/
# This process is divided into three basic step per app
# 1) Create an environment for the new app: conda create -n <app_environment_name>
# 2) Activate this environment: source activate <app_environment_name>
# 3) Install the app: conda install -c bioconda <app_name>
# 4) A list of packages was included in  the file condaapps.packs
# IMPORTANTE: Por um app em cada ambiente, ou seja um app por linha
echo "Installing Bioconda packages..."
source activate base # base environment
conda install h5py #  lets you store and manipulate huge amounts of numerical data from NumPy
while IFS= read -r line; do
  conda create -n "$line" "$line" -c bioconda -yq
  conda activate "$line" && "$line"
done < ~/repos/bioinfo/bioinfo_condaapps.packs

# Fim!
exit 0

# Configurar miniconda para outros usuários
#https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/admin-multi-user-install.html
#https://docs.anaconda.com/anaconda/install/multi-user/
#sudo groupadd conda
#sudo chgrp -R conda /home/kalabric/miniconda3
#sudo chmod 770 -R /home/kalabric/miniconda3
#sudo adduser bioinfo conda

# Instalar Guppy (ONT) - este procedimento não funcionou corretamente em 09/28/2023. Possívelmente, o guppy foi descontinuado e substituido por dorado!
#echo "Installing Guppy (ONT)..."
#sudo apt-get update
#sudo apt-get install wget lsb-release
#export PLATFORM=$(lsb_release -cs)
#wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
#echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
#sudo apt-get update
# 1. To install the .deb for Guppy, use the following command:
#sudo apt update
#sudo apt install ont-guppy
# This will install the GPU version of Guppy.
# or:
#sudo apt update
#sudo apt install ont-guppy-cpu
# To install the CPU-only version of Guppy

# Instalar Dorado (ONT)
# Link: https://github.com/nanoporetech/dorado
#echo "Installing Dorado (ONT)..."
#cd Downloads\
# 1) First method:
#wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz
#tar -xf dorado-0.3.4-linux-x64.tar.gz
# Copy all bin and lib file to your path
# 2) Second method:
# Requeriment:
sudo apt install cmake
# Clone and build
#cd repos
#git clone https://github.com/nanoporetech/dorado.git dorado
#cd dorado
#cmake -S . -B cmake-build
#cmake --build cmake-build --config Release -j
#ctest --test-dir cmake-build
#cmake --install cmake-build --prefix ${HOME}/bin

# Instalar DIAMOND - DIAMOND is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data
# downloading the tool
# wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz
# tar xzf diamond-linux64.tar.gz
# creating a diamond-formatted database file
# ./diamond makedb --in reference.fasta -d reference
# running a search in blastp mode
# ./diamond blastp -d reference -q queries.fasta -o matches.tsv
# running a search in blastx mode
# ./diamond blastx -d reference -q reads.fasta -o matches.tsv
# downloading and using a BLAST database
# update_blastdb.pl --decompress --blastdb_version 5 swissprot
# ./diamond prepdb -d swissprot
# ./diamond blastp -d swissprot -q queries.fasta -o matches.tsv


# Instalar Star
#https://github.com/alexdobin/STAR
# Get latest STAR source from releases
#wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz -P ~/Downloads
#tar -xzf ~/Downloads/2.7.9a.tar.gz -C ~/Downloads
# Compile
#cd ~/Downloads/STAR-2.7.9a/source/ # Reque cd na linha de comando
#make STAR
#cp STAR ~/bin/

# Instalar phred-phrap-consed
#http://www.phrap.org/consed/consed.html#howToGet
#cd ~/Downloads/
#mkdir phred-phrap-cap3-consed
#tar xzvf phred-dist-020425.c-acd.tar.Z --one-top-level
#cd phred-dist-020425.c-acd
#make
#cd ..
#find phred-dist-020425.c-acd/ -perm /a+x -exec cp {} ~/bin/ \;
#tar -xvf distrib.tar.Z --one-top-level
#cd distrib
#make
#cd ..
#find distrib/ -perm /a+x -exec cp {} ~/bin/ \;
#tar xzvf consed_linux.tar.gz --one-top-level
#cd consed_linux
#make
#cp phred ~/bin
#cp phredpar.dat ~/bin
#cd ~/Downloads/


#make
#cp phred ~/bin/phred-phrap-consed/bin

#nano .profile
	# set PATH and CONSED_HOME
	#CONSED_HOME=/home/kalabric/bin/phred-phrap-consed/bin
	#export PHRED_PARAMETER_FILE="$CONSED_HOME/phredpar.dat"
	#if [ -d "$CONSED_HOME" ] ; then
	#    PATH="$CONSED_HOME:$PATH"
	#fi
#make daev
#cp daev ~/bin/phred-phrap-consed/bin
#mkdir phd2fasta
#cp phd2fasta-acd-dist.130911.tar.gz phd2fasta
#cd phd2fasta
#tar xzvf phd2fasta-acd-dist.130911.tar.gz
#make
#cp phd2fasta ~/bin/phred-phrap-consed/bin
#mkdir phrap
#cp distrib.tar.Z phrap/
#cd phrap/
#tar xzvf distrib.tar.Z
# Compilar phrad-swat-cross_match
#make
#cp phrap ~/bin/phred-phrap-consed/bin
#cp cluster ~/bin/phred-phrap-consed/bin
#cp loco ~/bin/phred-phrap-consed/bin
#cp cross_match ~/bin/phred-phrap-consed/bin
#cp swat ~/bin/phred-phrap-consed/bin
#cp phrapview ~/bin/phred-phrap-consed/bin



# Instalar RepeatMasker
# Link: https://www.repeatmasker.org/RepeatMasker/
# Requesitos: python3, biblioteca pyton h5py, Cross_Match, RMBlast, TRF, hmmer (bioconda)
# wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz -P ~/Downloads
# wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz -P ~/Downloads
# wget https://www.dfam.org/releases/Dfam_3.2/families/Dfam.h5.gz -P ~/Downloads
# tar xzvf RepeatMasker-4.1.2-p1.tar.gz

# Instalar Glimmer
# Link: http://ccb.jhu.edu/software/glimmer/index.shtml
# wget http://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz -P ~/Downloads
# cd ~/Downloads
# tar xzvf glimmer302b.tar.gz
# cd glimmer3.02/src
# make

# Instalar Augustus
# Link: http://bioinf.uni-greifswald.de/augustus/
# sudo apt install augustus augustus-data augustus-doc

# Instalar gmap
# Link: https://github.com/juliangehring/GMAP-GSNAP
# wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2021-12-17.tar.gz
# tar -xvf gmap-gsnap-2021-12-17.tar.gz
# nano config.site
# .configure
# make
# make check
# make install
# gmap_build -d GRCh38_p13 -g data/GRCh38/GRCh38.p13.genome.fa.gz
### GMAP SAMTOOLS Pipeline
# gmapl -d GRCh38_p13 -A barcode01.fasta - t ${THREADS} -O ### Para genomas grandes usar gmapl
# gsnap -d GRCh38_p13 barcode01.fastq ou cat barcode01.fastq | gsnap -d GRCh38_p13
# samtools view -f 4 file.bam > unmapped.sam

# Instalar musket
# Link: https://musket.sourceforge.net/homepage.htm
# cd ~/Downloads
# wget https://sourceforge.net/projects/musket/files/latest/download/musket-1.1.tar -P ~/Downloads
# tar xzvf musket-1.1.tar
# cd musket-1.1
# make
# cp musket ~/bin

# Instalar flash
# Link: http://ccb.jhu.edu/software/FLASH/
# cd ~/Downloads
# wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11.tar.gz -P ~/Downloads
# tar xzvf FLASH-1.2.11.tar.gz
# cd FLASH-1.2.11
# make
# cp flash ~/bin

# Instalar velvet
# Link: https://github.com/dzerbino/velvet
# cd repos
# git clone https://github.com/dzerbino/velvet.git
# cd velvet/
# make 
# cp velvet* ~/bin
