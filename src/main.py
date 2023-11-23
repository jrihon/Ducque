import argparse
import sys
import time

from builder import builder             # Build the duplex
from transmute import transmute         # Convert a given pdb file to a custom json format
from randomise import randomise         # Output a random sequence
import process_CLI_inputs               # Exceptions and custom errors
import systemsDucque as SD              # Retrieves system information from the machine

def main():

    explanation = """
                                                          AUG   TCC        TCG   TTC        TAT   TCG
                             ____                         :::C A:::G      A:::G T:::A      A:::T C:::G
                            |  _ \\ _   _  ___ __ _ _   _  ___ G:::::A    G:::::T:::::G    G:::::G:::::C         :
                            | | | | | | |/ __/ _` | | | |/ _ \\ T:::::G  T:::::C A:::::T  G:::::A G:::::T   G   :G
                            | |_| | |_| | (_| (_| | |_| |  __/  A:::::TC:::::C   A:::::AT:::::A   A:::::AGA:  :A
                            |____/ \\__,_|\\___\\__, |\\__,_|\\___|   G::::CG::::G     G::::TA::::T     C::::TC::::T
                                                |_|               GCTC  AGTA       GGCA  CCTA       CCGA  TCTA

    Project to generate and customize DNA, RNA, XNA duplex molecules

    Designed and written by Doctorandus Rihon Jérôme.
    ______________________________________________________________________
    """

    t0 = time.time()

    # Check if the version of python is at least python3.8
    SD.version_checker()

    ## -------------------------------------------------- P A R S E  A R G U M E N T S --------------------------------------- ##
    options = argparse.ArgumentParser(
                                      description=explanation,
                                      add_help=False,
                                      formatter_class=argparse.RawTextHelpFormatter)

    options.add_argument("--build", type=argparse.FileType("r"),
            help="Build a nucleic acid duplex with a given sequence.")

    options.add_argument("--transmute", type=argparse.FileType("r"),
            help="Input a `*.pdb` of the nucleic acid you want to convert to `*.json`.\n")

    options.add_argument("--randomise", type=argparse.FileType("r"),
            help="Randomise a given sequence to-be-built.")

    options.add_argument("--xyz_pdb", type=argparse.FileType("r"),
            help="Convert an `*.xyz` file to a formatted `*.pdb` file.")

    options.add_argument("--gui", type=str, nargs='?', action="append", const="NO_FLAG",
            help="Call Ducque's GUI to access clickable actions (instead of using CLI).")

    options.add_argument("--help", action="help",
            help="Prompt Ducque's help message to appear")

    arguments = options.parse_args()

    if len(sys.argv) == 1:
        options.print_help()
        SD.exit_Ducque()

    else:
        # If we call the nucleic acid builder
        if arguments.build:
            NUCLEIC_ACID_LIST, COMPLEMENT, OUTFILE = process_CLI_inputs.build(arguments.build)

        # If we want to convert a pdb to a json file
        if arguments.transmute:
            PDB_FNAME, CHEMISTRY_T, MOIETY, DIHEDRALS, ANGLES, CONFORMATION, NUCLEOBASE = process_CLI_inputs.transmute(arguments.transmute)

        # If we want to call for a randomised sequence
        if arguments.randomise:
            CHEMISTRY_R, LENGTH_SEQUENCE, SEQUENCE, COMPL_SEQ, OUTFILE = process_CLI_inputs.randomise(arguments.randomise)

        # If we want to convert XYZ file to PDB
        if arguments.xyz_pdb:
            XYZ_FNAME, ATOM_ID, ATOMNAME_LIST = process_CLI_inputs.xyz_to_pdb(arguments.xyz_pdb)

        # If we want to use the GUI
        if arguments.gui :
            GUI_OPT = process_CLI_inputs.gui_module(arguments.gui)


## -------------------------------------------------- M A I N -------------------------------------------------- ##

    # Build nucleic acid duplex
    if arguments.build:
        SD.print_build()
        builder.Architecture(NUCLEIC_ACID_LIST, COMPLEMENT, OUTFILE)
        SD.print_time(time.time(), t0)
        sys.exit(0)

    # Convert pdb to json
    if arguments.transmute:
        SD.print_transmute(PDB_FNAME)
        transmute.Transmutation(PDB_FNAME, CHEMISTRY_T, MOIETY, DIHEDRALS, ANGLES, CONFORMATION, NUCLEOBASE)
        sys.exit(0)

    # Convert an xyz coordinate file to a pdb file
    if arguments.xyz_pdb:
        SD.print_xyz(XYZ_FNAME)
        transmute.convert_XYZ_to_PDB(XYZ_FNAME, ATOM_ID, ATOMNAME_LIST)
        sys.exit(0)

    # Output a randomised sequence
    if arguments.randomise:
        SD.print_rand()
        randomise.randomiser(CHEMISTRY_R, LENGTH_SEQUENCE, SEQUENCE, COMPL_SEQ, OUTFILE)
        sys.exit(0)

    # Use the GUI module. This gets loaded at the back because GUI users do not mind efficiency
    if arguments.gui :
        from dgui import ducqueGUI 
        print(" <^)    Ducque")
        print(" ( 3 )   [GUI]")
        ducqueGUI.select_window(GUI_OPT)




## ---------------------------------------- S T A R T   T H E   S H O W ---------------------------------------- ##
if __name__ == "__main__":
    main()
