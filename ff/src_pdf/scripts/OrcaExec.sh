#!/bin/bash


# Two ways of using this shell script : 
# $ ./OrcaExec.sh BASENAME
#
# $ bash OrcaExec.sh BASENAME



# You can call this and run this with taskspooler, a shell tool that runs jobs in the background (or simultaneously)
# $ tsp bash OrcaExec.sh BASENAME

# Make sure that you have changed to name of the pdb file you prompt into ORCA in the `template_Orca.inp`
# If necessary, also change the charge of the molecule (and possibly even the multiplicity, but for organic molecules this is almost always 1)


##### ----- START SCRIPT

# This works for Orca 5.0.2
# If ORCA is not installed in your $HOME directory, change this path.

# First argument that comes after this script. Should be the basename of the *.inp you are prompting to ORCA
INP=$1

$HOME/Orca-5.0.2-openmpi4.1.1/orca $INP.inp "--use-hwthread-cpus" > $INP.out # v5.0.2

