# Username

# Link: https://linuxize.com/post/how-to-create-users-in-linux-using-the-useradd-command/

# Editar o arquivo /etc/default/useradd para configurar a variável SHELL=/bin/bash

USERNAME = $1
sudo useradd $USERNAME
sudo passwd $USERNAME
