#!/bin/bash

DucqueDirectory=$(pwd) # Make sure we are in the correct directory


# See if by some sheer quirk the run file is not there or it is empty
if [[ -f "$DucqueDirectory/bin/Ducque" ]] && [[ -s "$DucqueDirectory/bin/Ducque" ]]; then
    chmod +x "$DucqueDirectory/bin/Ducque" # give executable permissions to the Ducque software
    echo "Ducque was successfully set up."
else
    echo "Ducque's software was not found."
    echo "Check if you are in the correct directory in order to run the setup.sh script"; echo
    echo "No such file or directory : $DucqueDirectory/bin/Ducque"
fi
