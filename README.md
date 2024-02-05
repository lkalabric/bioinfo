# Project bioinfo
<b>Repository Goal</b><br>
Our main goal is to make different scripts available to help users install bioinformatics packages, discuss their usage by providing example data and actually analyze this data. Simply creating this repository will be a great opportunity for all of us to learn by doing, by programming Linux's own bash scripting language, by running the packages themselves, or by creating scripts to run them individually or in a pipeline!

<b>Bio-Linux from scratch - Creating your own Bio-Linux Distro</b><br>
Since Bio-Linux (DistroWatch https://distrowatch.com/table.php?distribution=biolinux) was discontinued in 2015 (latest version 8.0.5), nobody has organized any other Linux distribution to replace it. I believe the distribution was discontinued because the Bioinfo data became so complex that we needed more powerful machines to do most of the analyses, good laptops became insufficient and CPUs gave way to GPUs. Even using a distribution focused on bioinformatics, considered one of the most complete, powerful, configurable, and easy-to-use systems, not all users would be able to run more complete pipelines. For those like me who are still interested in using a shared resource and learning from it, I hope this repository will also help you set up a Bio-Linux-like distro from scratch. Anyone interested in contributing is welcome!

<b>Required Operational System (OS)</b><br>
Any Linux distro can be used for this purpose, however, we strongly recommend all users to update and upgrade their Linux software before installing any package or running any of our scripts. All our experience was using Ubuntu Linux 22.04.2 LTS (Jammy Jellyfish). If you do not have sudo privileges, ask a sudo user to update/upgrade your system from time to time. 

<b>Local server</b><br>
To access our local server, we recommend you use Putty (available at https://putty.org/). Our IP address is <b>10.71.17.54</b>. The username and password are <b>bioinfo</b>!

To find the name of the Linux operating system and the kernel version you are running, type:<br>
$ lsb_release -a<br>

To get the updated software list for Ubuntu, enter:<br>
$ sudo apt-get update<br>
or<br>
$ sudo apt update<br>

To see available updates, run:<br>
$ sudo apt list --upgradable<br>

To update software(s) i.e. apply updates and patches on Ubuntu Linux, type the following:<br>
$ sudo apt-get upgrade<br>
or<br>
$ sudo apt upgrade<br>

<b>Installation</b><br>
Set up your filesystem, directory tree, PATH and start installing programming languages and bioinformatics packages. We intend to create different scripts for each step of the installation process.<br>
Important: Read the documentation available in each script before running it on your system.<br> 

First, let´s check if git package is installed and install it if not?<br>
$ dpkg -s git<br>
Note: git package should be installed in Ubuntu Linux 22.04.2 LTS (Jammy Jellyfish) by default. 

If this is not the case, run:<br>
$ sudo apt install git<br>

Second, let´s prepare our environment by cloning bioinfo repository from github into our local machine and be able to execute all scripts.<br>
$ git clone https://github.com/lkalabric/bioinfo.git repos/bioinfo    # clone bioinfo repo from github into the directory repos/<br>

Third, this is optional. We recommend the following directory tree to organize your files and directories:<br>
$ bash repos/bioinfo/create_directories.sh                            # create standardized directories and prepare your Linux filesystem

Suggestion of a directory tree<br>
bin/ - binary executable files<br>
data/ - your data files<br>
examples/ - example files<br>
logs/ - log files<br>
repos/ - github repositories<br>
results/ - your results<br>
scripts/ - commands created in different scripting languages<br>

Forth, let´s keep our repos/bioinfo and scripts up-to-date<br>
$ bash repos/bioinfo/git_scripts.sh bioinfo                           # git pull and copy all up-to-date scripts from repos/bioinfo to scripts/

Fifth, last but not least, add scripts/ to PATH<br> 
Note: The most common directories that hold executable programs are /bin, /sbin, /usr/sbin, /usr/local/bin and /usr/local/sbin<br>

To add scripts/ to PATH, run:<br>
$ export PATH="$HOME/scripts:$PATH"                                     # this exports $PATH<br>

To make the change permanent, you need to define the $PATH variable in a shell configuration file like ~/.bashrc.<br>
$ echo 'export PATH="$HOME/scripts:$PATH" # add scripts/ to PATH' >> ~/.bashrc  # appends the export to ~/.bashrc file<br>

After saving the file, run the export take effect:<br>
$ source ~/.bashrc

To test if everything is working so far, let´s do this:
$ git_scripts.sh bioinfo                                             # This command pull all files from bioinfo repo, copy all *.sh files to the script/ diretory and changes their mode to executable

Most used packages from different repository
- Perl (CPAN)
- Python3
$ sudo apt-get install python3
- Conda or Miniconda
Link: https://docs.conda.io/en/latest/miniconda.html
$ mkdir -p ~/miniconda3
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o ~/miniconda3/miniconda.sh
$ bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
$ rm -rf ~/miniconda3/miniconda.sh
$ ~/miniconda3/bin/conda init bash
$ ~/miniconda3/bin/conda init zsh
- R (BioConductor)
$ sudo apt install r-cran-littler

To make it simpler, one can install all packages runnig:<br>
$ install_linuxpacks.sh -i *.packs<br>
$ install_thirdparty.sh

<b>Pipeline ready-to-use steps (MinIon & Illumina)</b><br>
- organize_filesystem.sh
- basecalling_demux.sh (MinIon only)
- quality_control.sh
- quality_filter.sh
- assembly_de_novo.sh <-spades / -velvet><br>
Note: We intend to develop a code to generate full pipelines using different software and take advantage of the codes already written so far.

<b>Database</b><br>
- Viral Genomes > Complete RefSeq release of viral sequences (available at https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) 
