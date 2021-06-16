import argparse
import sys
import os
from time import time

import transmute    # Convert a given pdb file to a custom json format
import labyrinth    # Build the duplex
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

options.add_argument("--randomise",
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
            print("One or more of the nucleotides in the given sequence is invalid. Please check your input file.\n\n\n")
            options.print_help()
            sys.exit(0)


    # If we want to convert a pdb to a json file
    if arguments.transmute:
        # Read the input file and create a list object of the sequence. Removes any whitespace.
        fileTransmute = list(map(lambda x: x.strip(), arguments.transmute.readlines()))

        # Check if amount of inputs are valid
        list_of_valid_flags = ["--pdb", "--id"]
        if len(fileTransmute) != len(list_of_valid_flags):
            print("Only two arguments are required, please check our input file.\n\n\n")
            options.print_help()
            sys.exit(0)

        # Check if flags are valid
        for argument in fileTransmute:
            arg = argument.split()
            # If a given flag is not a valid one, shut it down
            if not arg[0] in list_of_valid_flags:
                print(arg[0])
                print("\n\nOne or more prompted flags are not valid, please check our input file.\n\n\n")
                options.print_help()
                sys.exit(0)

            #If the argument is a numerical value, shut it down
            if arg[1].isdigit():
                print("\n\nThe following argument has been passed: ", arg[1])
                print("Arguments can not be of a numerical value, please check our input file.\n\n\n")
                options.print_help()
                sys.exit(0)

            # Save the name of the pdb file as a variable
            if arg[0] == "--pdb":
                pdb_file = arg[0]

            # Save the identifier of the nucleic acid as a variable
            if arg[0] == "--id":
                identifier = arg[1]


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
        print('Reading nucleic acid molecule ...')
        transmute.Transmutation(pdb_file, identifier)


    # Build nucleic acid duplex
    if arguments.Daedalus:
        labyrinth.Architecture(fileDaedalus)

    print('Time spent: %.5f seconds.' % (time() - t0))




## ---------------------------------------- S T A R T   T H E   S H O W ---------------------------------------- ##
if __name__ == '__main__':
    main()






## ----------------------------------- FUNCTIONS THAT ARE NOT IN USE ANYMORE ----------------------------------- ##
#    # Convert json to pdb
#    if arguments.json and transmute.json_is_linker(arguments.json):
#        print('Reading linker .json ...')
#        transmute.linker_to_pdb(arguments.json)
#        sys.exit()
#    elif arguments.json:
#        print('Reading nucleic acid .json ...')
#        transmute.NA_to_pdb(arguments.json)
#        sys.exit(0)
