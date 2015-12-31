#!/usr/bin/env bash

#Check if virtualenv is installed and install
pip install virtualenv

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

#Re-source bash_profile
source ~/.bash_profile

#Deactivate virtualenv
deactivate

#Backup to original directory
cd ..
cd ..

