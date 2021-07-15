import argparse
import sys
import os
import time

import transmute    # Convert a given pdb file to a custom json format
import labyrinth    # Build the duplex
import randomise    # Output a random sequence
import fundaments   # Exceptions or custom errors

explanation = """
                                                               AUG   TCC        TCG   TTC        TAT   TCG
                         ____                 _       _        :::C A:::G      A:::G T:::A      A:::T C:::G
                        |  _ \  __ _  ___  __| | __ _| |_   _ ___  G:::::A    G:::::T:::::G    G:::::G:::::C         :
                        | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|  T:::::G  T:::::C A:::::T  G:::::A G:::::T   G   :G
                        | |_| | (_| |  __/ (_| | (_| | | |_| \__ \   A:::::TC:::::C   A:::::AT:::::A   A:::::AGA:  :A
                        |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/    G::::CG::::G     G::::TA::::T     C::::TC::::T
                                                                       GCTC  AGTA       GGCA  CCTA       CCGA  TCTA

Project to generate and customize DNA, RNA, XNA duplex molecules

Designed and written by Doctorandus Rihon Jérôme.\n
______________________________________________________________________
"""

t0 = time.time()

## -------------------------------------------------- P A R S E  A R G U M E N T S --------------------------------------- ##
options = argparse.ArgumentParser(description=explanation,
                                  add_help=False,
                                  formatter_class=argparse.RawTextHelpFormatter)

options.add_argument("--transmute", type=argparse.FileType("r"),
        help="Input the pdb of the nucleic acid you want to convert to json.\n")

options.add_argument("--randomise", type=argparse.FileType("r"),
        help="Given a set of parameters, write out a random sequence that can be prompted to --Daedalus.")

options.add_argument("--Daedalus", type=argparse.FileType("r"),
        help="Ask the software to build a nucleic acid duplex with a given sequence.")

options.add_argument("--help", action="help",
        help="Prompt Daedalus's help message to appear")

arguments = options.parse_args()


if len(sys.argv) == 1:
    options.print_help()
    sys.exit(0)

else:
    # If we call the nucleic acid builder
    if arguments.Daedalus:
        nucleic_acid_list, complement = fundaments.daedalus(arguments.Daedalus, options)

    # If we want to convert a pdb to a json file
    if arguments.transmute:
        pdb_file, identifier, moiety, dihedrals, angles, conformation = fundaments.transmute(arguments.transmute, options)

    # If we want to call for a randomised sequence
    if arguments.randomise:
        chemistry, length_sequence, sequence = fundaments.randomise(arguments.randomise, options)


# Try and see if any of the options are used together. 
# Daedalus will not allow this to happen for the reason that I don't feel like complicating stuff too much.
try:
    if arguments.Daedalus and arguments.transmute:
        raise fundaments.InputExclusivity

except fundaments.InputExclusivity:
    print("InputExclusivity: These flags are mutually exclusive; --transmute --Daedalus ")
    sys.exit(0)


## -------------------------------------------------- M A I N -------------------------------------------------- ##
def main():

    # Print the Daedalus prompt
    print(explanation)

    # Convert pdb to json
    if arguments.transmute:
        transmute.Transmutation(pdb_file, identifier, moiety, dihedrals, angles, conformation)
        print("Converting " + pdb_file + " to a json file.")

    # Build nucleic acid duplex
    if arguments.Daedalus:
        print("Daedalus - Nucleic Acid Builder initiated! Building sequence ...\n")
        labyrinth.Architecture(nucleic_acid_list, complement)

    # Output a randomised sequence
    if arguments.randomise:
        print("Randioli randioli, what is the spaghetolli?")
        randomise.randomiser(chemistry, length_sequence, sequence)

    print("                         Time spent: %.5f seconds." % (time.time() - t0))


## ---------------------------------------- S T A R T   T H E   S H O W ---------------------------------------- ##
if __name__ == "__main__":
    main()
