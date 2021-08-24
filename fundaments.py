import labyrinth_func_tools3 as LFT3
import sys

""" When the user prompts the wrong values or flags, this python script intercepts most errors that happen at the start """

class InputExclusivity(Exception):
    Except = " These flags are mutually exclusive; --transmute     --Daedalus "



def check_if_nucleotides_are_valid(input_sequence : list) -> bool:
    """ Check if any of the prompted nucleotides is not valid. """
    # Retrieve the keys of the dictionary from which we parse the data
    keys_of_dict = LFT3.codex_acidum_nucleicum.keys()
    # Check if any of the prompted nucleotides is not found in the sequence

    for NA in input_sequence:
        if NA not in keys_of_dict:
            print("One or more of the nucleotides in the given sequence is invalid. Please check your input file : " + NA + "\n\n")
            return False

    return True


def daedalus(DaedalusInput, options):

    list_of_valid_flags = ["--sequence", "--complement"]
    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileDaedalus = list(map(lambda x: x.strip(), DaedalusInput.readlines()))

    # Check if amount of inputs are valid
    if len(fileDaedalus) != len(list_of_valid_flags):
        print("Only two arguments are required, please check our input file.\n\n\n")
        options.print_help()
        sys.exit(0)

    for argument in fileDaedalus:
        arg = argument.split()

        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print("\n\nThe following flag is invalid : " + arg[0] + ". Please check your input file.\n\n\n")
            #options.print_help()
            sys.exit(0)

        if arg[0] == "--sequence":
            nucleic_acid_list = list(map(lambda x: x.strip(","), arg[1:]))

        if arg[0] == "--complement":
            # If there is a input possibility at index 2, meaning more than one string have been inputted, then get the entire string as a list variable.
            # If there is not an input possibility at index 2, this means there is only one input after the flag available and that means it is just a string.
            try:
                arg[1]
            except:
                print("You have forgotten to include an argument for the --complement flag! Please add this to your input file.")
                sys.exit(0)

            try:
                arg[2]
            except:
                complement = arg[1]
            else:
                complement = arg[1:]

    # If a given nucleotide is not in the list of valid nucleotides, stop the program
    if not check_if_nucleotides_are_valid(nucleic_acid_list):
        #options.print_help()
        sys.exit(0)

    return nucleic_acid_list, complement


def transmute(TransmuteInput, options):

    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileTransmute = list(map(lambda x: x.strip(), TransmuteInput.readlines()))

    # Check if amount of inputs are valid
    list_of_valid_flags = ["--pdb", "--id", "--moiety", "--conformation", "--dihedrals", "--bondangles"]
    if not len(fileTransmute) == len(list_of_valid_flags) and not len(fileTransmute) == (len(list_of_valid_flags) - 1):
        print("Only five/six arguments are required, please check your input file.\n\n\n")

        options.print_help()
        sys.exit(0)

    for argument in fileTransmute:
        arg = argument.split()
        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print("\n\nThe following flag is invalid : " + arg[0] + ". Please check your input file.\n\n\n")
            #options.print_help()
            sys.exit(0)

        # Save the prompted arguments according to the respective flags
        if arg[0] == "--pdb":
            pdb_file = arg[1]

        if arg[0] == "--id":
            identifier = arg[1]

        if arg[0] == "--conformation":
            conformation = arg[1]

        if arg[0] == "--moiety":
            moiety = arg[1]

        if arg[0] == "--dihedrals":
            dihedrals = arg[1:]

        if arg[0] == "--bondangles":
            angles = arg[1:]

    return pdb_file, identifier, conformation, moiety, dihedrals, angles


def randomise(RandomiseInput, options):

    fileRandomise = list(map(lambda x: x.strip(), RandomiseInput.readlines()))

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
            #options.print_help()
            sys.exit(0)

        # Save the prompted arguments according to their respective flags
        if arg[0] == "--chemistry":
            # If multiple chemistries are prompted
            if len(arg) > 2:
                chemistry = arg[1:]
            # If one chemistry is prompted
            elif len(arg) == 2:
                chemistry = arg[1]
        if arg[0] == "--length":
            length_sequence = int(arg[1])
            sequence = None

        if arg[0] == "--sequence":
            sequence = arg[1:]
            length_sequence = 0

    # If the variable chemistry has not been defined, then
    try:
        chemistry
    except NameError:
        print("No chemistry type was prompted! Revise your input file. \n")
        #options.print_help()
        sys.exit(0)

    return chemistry, length_sequence, sequence


def xyz_to_pdb(ConversionInput, options):

    fileConversion = list(map(lambda x: x.strip(), ConversionInput.readlines()))

    # Check list of valid inputs 
    list_of_valid_flags = ["--xyz", "--atomID", "--atomname_list"]
    if len(fileRandomise) != 3:
        print("Only three arguments are required at one time. Please check your input file. \n\n\n ")
        options.print_help()
        sys.exit(0)

    # Check if flags are valid
    for argument in ConversionInput:
        arg = argument.split()
        # If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print("\n\nThe following flag is invalid : " + arg[0] + ". Please check your input file.\n\n\n")
            #options.print_help()
            sys.exit(0)

        # Save the prompted arguments according to their respective flags
        if arg[0] == "--xyz":
            fname_xyz = arg[1]
            try:
                os.path.isfile(fname_xyz)
            except FileNotFoundError:
                print("File " + fname_xyz + "was not found. Please revise either its name or its location.\n")
                sys.exit(0)

        if arg[0] == "--atomID":
            atomID = arg[1]
            if len(atomID) != 3:
                print("The argument '--atomID' requires three characters in its name.\n"
                        "Technically it does not, but I made Daedalus require it, so deal with it.\n")

        if arg[0] == "--atomname_list":
            atomname_list = arg[1:]

    return fname_xyz, atomID, atomname_list

