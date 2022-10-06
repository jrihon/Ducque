import argparse
import sys
import time

import labyrinth    # Build the duplex
import transmute    # Convert a given pdb file to a custom json format
import randomise    # Output a random sequence
import pre_processing   # Exceptions and custom errors
import sysDaedalus  # Retrieves system information from the machine

def main():
    explanation = """
                                                                   AUG   TCC        TCG   TTC        TAT   TCG
                             ____                 _       _        :::C A:::G      A:::G T:::A      A:::T C:::G
                            |  _ \  __ _  ___  __| | __ _| |_   _ ___  G:::::A    G:::::T:::::G    G:::::G:::::C         :
                            | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|  T:::::G  T:::::C A:::::T  G:::::A G:::::T   G   :G
                            | |_| | (_| |  __/ (_| | (_| | | |_| \__ \   A:::::TC:::::C   A:::::AT:::::A   A:::::AGA:  :A
                            |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/    G::::CG::::G     G::::TA::::T     C::::TC::::T
                                                                           GCTC  AGTA       GGCA  CCTA       CCGA  TCTA

    Project to generate and customize DNA, RNA, XNA duplex molecules

    Designed and written by Doctorandus Rihon Jérôme.
    ______________________________________________________________________
    """

    t0 = time.time()

    # Check if the version of python is at least python3.6
    sysDaedalus.version_checker()

    ## -------------------------------------------------- P A R S E  A R G U M E N T S --------------------------------------- ##
    options = argparse.ArgumentParser(description=explanation,
                                      add_help=False,
                                      formatter_class=argparse.RawTextHelpFormatter)

    options.add_argument("--Daedalus", type=argparse.FileType("r"),
            help="Ask the software to build a nucleic acid duplex with a given sequence.")

    options.add_argument("--transmute", type=argparse.FileType("r"),
            help="Input the *.pdb of the nucleic acid you want to convert to *.json.\n")

    options.add_argument("--randomise", type=argparse.FileType("r"),
            help="Given a set of parameters, write out a random sequence that can be prompted to --Daedalus.")

    options.add_argument("--xyz_pdb", type=argparse.FileType("r"),
            help="Convert the inputted *.xyz file to a properly formatted *.pdb file.")

    options.add_argument("--help", action="help",
            help="Prompt Daedalus's help message to appear")

    arguments = options.parse_args()


    if len(sys.argv) == 1:
        options.print_help()
        sysDaedalus.exit_Daedalus()

    else:
        # If we call the nucleic acid builder
        if arguments.Daedalus:
            NUCLEIC_ACID_LIST, COMPLEMENT, OUTFILE = pre_processing.daedalus(arguments.Daedalus, options)

        # If we want to convert a pdb to a json file
        if arguments.transmute:
            PDB_FNAME, CHEMISTRY_T, MOIETY, DIHEDRALS, ANGLES, CONFORMATION = pre_processing.transmute(arguments.transmute, options)

        # If we want to call for a randomised sequence
        if arguments.randomise:
            CHEMISTRY_R, LENGTH_SEQUENCE, SEQUENCE, COMPL_SEQ, OUTFILE = pre_processing.randomise(arguments.randomise, options)

        # If we want to convert XYZ file to PDB
        if arguments.xyz_pdb:
            XYZ_FNAME, ATOM_ID, ATOMNAME_LIST = pre_processing.xyz_to_pdb(arguments.xyz_pdb, options)


    # Try and see if any of the options are used together. 
    # Daedalus will not allow this to happen for the reason that I don't feel like complicating stuff too much.
    try:
        if arguments.Daedalus and arguments.transmute:
            raise pre_processing.InputExclusivity

    except pre_processing.InputExclusivity:
        print("InputExclusivity: These flags are mutually exclusive; --transmute --Daedalus ")
        sysDaedalus.exit_Daedalus()


    ## -------------------------------------------------- M A I N -------------------------------------------------- ##


    # Print the Daedalus prompt
#    print(explanation)

    # Build nucleic acid duplex
    if arguments.Daedalus:
        labyrinth.Architecture(NUCLEIC_ACID_LIST, COMPLEMENT, OUTFILE)

        print(f"                         Time spent: %.5f seconds." % (time.time() - t0))

    # Convert pdb to json
    if arguments.transmute:
        print("-----------------------------------------------------------")
        print(f"Converting {PDB_FNAME} to a json file.")
        transmute.Transmutation(PDB_FNAME, CHEMISTRY_T, MOIETY, DIHEDRALS, ANGLES, CONFORMATION)

    # Convert an ORCA xyz-formatted molecule file to a pdb file
    if arguments.xyz_pdb:
        print("-----------------------------------------------------------")
        print(f"Converting {XYZ_FNAME} to a .pdb format file.\n")
        transmute.convert_XYZ_to_PDB(XYZ_FNAME, ATOM_ID, ATOMNAME_LIST)

    # Output a randomised sequence
    if arguments.randomise:
        print("-----------------------------------------------------------")
        print("Randomisation of the given inputs!\n")
        randomise.randomiser(CHEMISTRY_R, LENGTH_SEQUENCE, SEQUENCE, COMPL_SEQ, OUTFILE)



## ---------------------------------------- S T A R T   T H E   S H O W ---------------------------------------- ##
if __name__ == "__main__":
    main()
