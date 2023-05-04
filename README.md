# biolinux
<b>Biolinux from scratch - Creating your own Bio-Linux Distro</b><br>
Since Bio-Linux (DistroWatch https://distrowatch.com/table.php?distribution=biolinux) was discontinued in 2015 (latest version 8.0.5), nobody has organized any other Linux distribution to replace it. I believe the distribution was discontinued because the Bioinfo data became so complex that we needed more powerful machines to do most of the analyses, good laptops became insufficient and CPUs gave way to GPUs. Even using a distribution focused on bioinformatics, considered one of the most complete, powerful, configurable and easy-to-use systems, not all users would be able to run more complete pipelines. For those like me who are still interested in using such a shared resource and learning by doing, I hope there is still hope. Anyone interested to contribuite is welcome!

<b>Repository Goal</b><br>
Our main goal is to make different scripts available to help users install bioinformatics packages, discuss their usage by providing example data and actually analyze this data. Simply creating this repository will be a great opportunity for all of us to learn by doing, by programming Linux's own bash scripting language, by running the packages themselves, or by creating scripts to run them individually or in a pipeline!

<b>Required Operational System</b><br>
Any Linux distro can be used for this purpose, however, we strongly recommend all users to update and update their Linux software before installing any packages or running any of our scripts. All our experience was using Ubuntu Linux 22.04.2 LTS (Jammy Jellyfish). 

To find the name of the Linux operating system and the kernel version you are running, type:
$ lsb_release -a

Get updated software list for Ubuntu, enter:
$ sudo apt-get update
or
$ sudo apt update

To see available updates, run:
$ sudo ap list --upgradable

Update software(s) i.e. apply updates and patches on Ubuntu Linux, type the following:
$ sudo apt-get upgrade
or
$ sudo apt upgrade

<b>Installation</b><br>
Programming languages, package repositories and bioinformatics packages themselves will be installed by different scripts. Read their documentation before running them on your system. Most of the time, we will only install packages that you don't already have installed on your system and that are necessary for our bioinformatics training!

First, let´s check if git package is installed and install it if not?
$ dpkg -s git
NOTE: git package ins installed in Ubuntu Linux 22.04.2 LTS (Jammy Jellyfish) by default. If this is not the case, type:
$ sudo apt install git

Second, let´s prepare our enviroment to clone biolinux git repository into our local machine and be able to execute all scripts. 
$ mkdir scripts                                           # create a directory named scripts/
$ mkdir repos                                             # create a directory named repos/
$ cd repos                                                # change to directory repos
$ git clone https://github.com/lkalabric/biolinux.git     # clone biolinux git repo into the directory repos/biolinux/
$ cp biolinux/*.sh ~/scripts/                             # copy all scripts to the directory scripts/
$ cd                                                      # change to your home directory
$ chmod +x scripts/*                                      # change scripts mode to executable

Third, this is optional. We recommend the following directory tree to organize your files and directories:
$ create_folders.sh
bin/ - 
repos/ - 
scripts/ - 
examples/ -
data/ -
results/ -


