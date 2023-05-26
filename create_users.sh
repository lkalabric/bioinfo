#!/bin/bash

# Username
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 26 MAY 2023
# objetive: Create users in Linux
# syntax: ./create_users.sh <username>
# link: https://linuxize.com/post/how-to-create-users-in-linux-using-the-useradd-command/

# Edit file sudo nano /etc/default/useradd to SHELL=/bin/bash

# Receive the username from the command line and store in a variable
USERNAME = $1

# Basic Linux commands
sudo useradd $USERNAME
sudo passwd $USERNAME
