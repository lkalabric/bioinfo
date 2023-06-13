# Atualizar o python 3.9
# https://www.pythoncentral.io/how-to-update-python/
# https://gist.github.com/patrickmmartin/5b6b2ddecd29d6a1b2ffee2d8eea31ec

# Instalar e configurar o Miniconda
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -P ~/Downloads
# bash ~/Downloads/Miniconda3-latest-Linux-x86_64.sh
# Instalar pacotes Bioconda
# https://anaconda.org/
# source activate base
# conda install h5py
# conda install -c bioconda iqtree
# conda install -c bioconda trimmomatic
# conda install -c bioconda cutadapt # Requer gcc 9
# conda install -c bioconda nanofilt
# conda install -c bioconda hmmer

# Configurar miniconda para outros usuários
#https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/admin-multi-user-install.html
#https://docs.anaconda.com/anaconda/install/multi-user/
#sudo groupadd conda
#sudo chgrp -R conda /home/kalabric/miniconda3
#sudo chmod 770 -R /home/kalabric/miniconda3
#sudo adduser bioinfo conda

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
#mkdir /home/kalabric/bin/phred-phrap-consed/
#cd ~/Downloads/phred-phrap-cap3-consed
#mkdir consed
#cp consed_linux.tar.gz consed
#cd consed/
#tar xzvf consed_linux.tar.gz
#./installConsed.perl consed_linux32bit /home/kalabric/bin/phred-phrap-consed/
#cd ~/Downloads/phred-phrap-cap3-consed
#mkdir phred
#cp phred-dist-020425.c-acd.tar.Z phred
#cd phred
#tar xzvf phred-dist-020425.c-acd.tar.Z
#Refer to the resulting file called 'INSTALL' for instructions on how to build and install phred
#make
#cp phred ~/bin/phred-phrap-consed/bin
#cp phredpar.dat ~/bin
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
