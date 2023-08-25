#!/bin/bash


# For example, the tarball (similar to zipfile) I downloaded was msms.i86_64Linux2.2.6.1
echo "Change the name of the LinuxX_X.X.X, on the BOTH LINES to the correct tar file name."
echo "Then remove or comment out the exit command from the script when succesful."
exit

# Make a new directory in your $HOME directory
mkdir $HOME/msms

# Extract the tarball file.
tar xvzf $HOME/Downloads/msms_i86_64LinuxX_X.X.X.tar.gz --directory $HOME/msms
cd $HOME/msms

# Copy the binary to a more manageable name, like `msms`
cp msms.x86_64LinuxX.X.X.X msms
