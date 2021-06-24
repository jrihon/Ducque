import argparse
import sys
import os
from time import time

import transmute    # Convert a given pdb file to a custom json format
import labyrinth    # Build the duplex
import randomise    # Output a random sequence
import fundaments   # Exceptions or custom errors

explanation = """
Project to generate and customize DNA, RNA, XNA duplex molecules

Designed and written by Doctorandus Rihon Jérôme.\n
______________________________________________________________________
"""
t0 = time()
## -------------------------------------------------- P A R S E  A R G U M E N T S --------------------------------------- ##
options = argparse.ArgumentParser(description=explanation,
                                  add_help=False,
                                  formatter_class=argparse.RawTextHelpFormatter)

## Let's not make this argument a required one, since we also might want to convert json to pdb
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
        # Read the input file and create a list object of the sequence. Removes any whitespace.
        fileDaedalus = list(map(lambda x: x.strip(), arguments.Daedalus.readline().strip("\n").split(",")))

        # If a given nucleotide is not in the list of valid nucleotides, stop the program
        if not fundaments.check_if_nucleotides_are_valid(fileDaedalus):
            options.print_help()
            sys.exit(0)


    # If we want to convert a pdb to a json file
    if arguments.transmute:
        # Read the input file and create a list object of the sequence. Removes any whitespace.
        fileTransmute = list(map(lambda x: x.strip(), arguments.transmute.readlines()))

        # Check if amount of inputs are valid
        list_of_valid_flags = ["--pdb", "--id", "--moiety", "--dihedrals", "--bondangles"]
        if len(fileTransmute) != len(list_of_valid_flags):
            print("Only five arguments are required, please check our input file.\n\n\n")
            options.print_help()
            sys.exit(0)

        # Check if flags are valid
        for argument in fileTransmute:
            arg = argument.split()
            # If a given flag is not a valid one, shut it down
            if not arg[0] in list_of_valid_flags:
                print("\n\nThe following flag is invalid : " + arg[0] + ". Please check your input file.\n\n\n")
                options.print_help()
                sys.exit(0)

            # Save the prompted arguments according to the respective flags
            if arg[0] == "--pdb":
                pdb_file = arg[1]

            if arg[0] == "--id":
                identifier = arg[1]

            if arg[0] == "--moiety":
                moiety = arg[1]

            if arg[0] == "--dihedrals":
                dihedrals = arg[1:]

            if arg[0] == "--bondangles":
                angles = arg[1:]

    # If we want to call for a randomised sequence
    if arguments.randomise:
        fileRandomise = list(map(lambda x: x.strip(), arguments.randomise.readlines()))

        # Check list of valid inputs 
        list_of_valid_flags = ["--chemistry", "--length", "--sequence"]
        if len(fileRandomise) != 2:
            print("Only two arguments are required at one time, please check our input file.\n\n\n")
            options.print_help()
            sys.exit(0)

        # Check if flags are valid
        for argument in fileRandomise:
            arg = argument.split()
            # If a given flag is not a valid one, shut it down
            if not arg[0] in list_of_valid_flags:
                print("\n\nThe following flag is invalid : " + arg[0] + ". Please check your input file.\n\n\n")
                options.print_help()
                sys.exit(0)

            # Save the prompted arguments according to their respective flags
            if arg[0] == "--chemistry":
                # If multiple chemistries are prompted
                if len(arg) > 2:
                    chemistry = arg[1:]
                # If one chemistry is prompted
                elif len(arg) == 2:
                    chemistry = arg[1]
            # If the variable chemistry has not been defined, then
            try:
                chemistry
            except NameError:
                print("No chemistry type was prompted! Revise your input file. \n")
                options.print_help()
                sys.exit(0)

            if arg[0] == "--length":
                length_seq = int(arg[1])
                sequence = None

            if arg[0] == "--sequence":
                sequence = arg[1:]
                length_seq = 0

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

    # Convert pdb to json
    if arguments.transmute:
        transmute.Transmutation(pdb_file, identifier, moiety, dihedrals, angles)
        print("Converting " + pdb_file + " to a json file.")

    # Build nucleic acid duplex
    if arguments.Daedalus:
        print("Bulding sequence ...\n")
        labyrinth.Architecture(fileDaedalus)

    # Output a randomised sequence
    if arguments.randomise:
        print("Randioli randioli, what is the spaghetolli?")
        randomise.randomiser(chemistry, length_seq, sequence)

    print("Time spent: %.5f seconds." % (time() - t0))

## ---------------------------------------- S T A R T   T H E   S H O W ---------------------------------------- ##
if __name__ == "__main__":
    main()
