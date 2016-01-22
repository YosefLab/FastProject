#!/usr/bin/env bash

echo
echo "Installing into: $(pwd)/FastProject"
echo

# Miniconda doesn't work for directory structures with spaces
if [[ $(pwd) == *" "* ]]
then
    echo "ERROR: Cannot install into a directory with a space in its path" >&2
    echo "Exiting..."
    echo
    exit 1
fi

# Test if new directory is empty.  Exit if it's not
if [ "$(ls -A $(pwd)/FastProject)" ]; then
    echo "ERROR: Directory is not empty" >&2
    echo "If you want to install into $(pwd)/FastProject, "
    echo "clear the directory first and run this script again."
    echo "Exiting..."
    echo
    exit 1
fi

curl "https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh" -o Miniconda_Install.sh

bash Miniconda_Install.sh -b -f -p FastProject

PATH="$(pwd)/FastProject/bin":$PATH

conda install numpy scipy matplotlib scikit-learn pandas -y

pip install FastProject-*.*.*.tar.gz

cd FastProject
mkdir Scripts
ln -s ../bin/fastproject Scripts/fastproject

echo "FastProject script installed to $(pwd)/Scripts"
echo
echo "Add folder to path by appending to .bash_profile?"
read -p "[y/n] >>> " -r
echo
if [[ ($REPLY == "yes") || ($REPLY == "Yes") || ($REPLY == "YES") ||
    ($REPLY == "y") || ($REPLY == "Y")]]
then
    echo "export PATH=\"$(pwd)/Scripts\":\$PATH" >> ~/.bash_profile
    echo "Your .bash_profile was updated."
    echo "Restart the terminal for the change to take effect"
else
    echo "Your PATH was not modified."
fi

# Cleanup
cd ..
rm Miniconda_Install.sh

echo
echo "FastProject successfully installed."


