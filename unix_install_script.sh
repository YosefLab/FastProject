#!/usr/bin/env bash
set -e

#Install virtualenv if not installed
pip install virtualenv --user

# For some reason, the user directory might not be on the path.
# Add it there anyways, just in case
export PATH=$PATH:$HOME/Library/Python/2.7/bin/

#Create virtualenv for FastProject
virtualenv FastProject
source FastProject/bin/activate
pip install FastProject-0.9.2.tar.gz

#Create a symlink to the script
cd FastProject
mkdir Scripts
cd Scripts
ln -s ../bin/fastproject fastproject

#Add to path
echo "" >> ~/.bash_profile
echo "# Added by FastProject Installer" >> ~/.bash_profile
echo "export PATH=\$PATH:$(pwd)" >> ~/.bash_profile

#Deactivate virtualenv
echo "Deactivate virtualenv"
deactivate

#Re-source bash_profile
echo "Re-loading profile"
source ~/.bash_profile

#Backup to original directory
cd ..
cd ..

