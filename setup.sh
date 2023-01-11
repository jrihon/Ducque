#!/bin/bash

DucqueDirectory=$(pwd) # Make sure we are in the correct directory


# See if by some sheer quirk the run file is not there
if [[ -f "$DucqueDirectory/bin/Ducque" ]] && [[ -s "$DucqueDirectory/bin/Ducque" ]]; then
    chmod +x "$DucqueDirectory/bin/Ducque" # give executable permissions to the Ducque software
fi
