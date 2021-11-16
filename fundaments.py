import labyrinth_func_tools3 as LFT3
import sys
import os

""" When the user prompts the wrong values or flags, this python script intercepts most errors that happen at the start """

def print_divide_between_command_and_output():
    """ This function exists solely to split the Daedalus call command and its output."""
    print("-----------------------------------------------------------")


class InputExclusivity(Exception):
    Except = " These flags are mutually exclusive; --transmute     --Daedalus "


def remove_trailing_whitespace(fileList : list) -> list:

    # intialise new list
    newList = []
    # Remove all trailing whitespace, cannot depend on only the last line being whitespace
    for i in fileList:
        # Truthy statement check. If this is not empty, append it to the new list
        if i.strip():
            newList.append(i)

    return newList

def check_if_nucleotides_are_valid(input_sequence : list) -> bool:
    """ Check if any of the prompted nucleotides is not valid. """
    # Retrieve the keys of the dictionary from which we parse the data
    keys_of_dict = LFT3.codex_acidum_nucleicum.keys()
    # Check if any of the prompted nucleotides is not found in the sequence

    for NA in input_sequence:
        if NA not in keys_of_dict:
            print_divide_between_command_and_output()
            print(f"One or more of the nucleotides in the given sequence is invalid. Please check your input file : {NA} \n\n")
            return False

    return True


def check_if_chemistry_is_valid(chemistry : str) -> list:
    """ If the chemistry is invalid, stop the script.
        If the chemistry is actually valid, return the keys of the dictionary as a list. """

    from transmute_func_tools import nucleoside_dict

    # Get the value from the transmute_func_tools nucleoside dictionary, to get the full name of the chemistry
    try:
        NUC_ID = nucleoside_dict[chemistry.upper()]
    except KeyError:
        print(f"The following key does not exist in the dictionary : {complement.upper()}.\nPlease revise your inputs.\n")
        sys.exit(1)

    return NUC_ID


def daedalus(DaedalusInput, options):

    list_of_valid_flags = ["--sequence", "--complement"]
    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileDaedalus = remove_trailing_whitespace(list(map(lambda x: x.strip(), DaedalusInput.readlines())))

    # Check if amount of inputs are valid
    if len(fileDaedalus) != len(list_of_valid_flags):
        print_divide_between_command_and_output()
        print("Only two arguments are required, please check our input file.\n\n\n")
        options.print_help()
        sys.exit(0)

    for argument in fileDaedalus:
        arg = argument.split()

        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
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
                print_divide_between_command_and_output()
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
    fileTransmute = remove_trailing_whitespace(list(map(lambda x: x.strip(), TransmuteInput.readlines())))

    # Check if amount of inputs are valid
    list_of_valid_flags = ["--pdb", "--id", "--moiety", "--conformation", "--dihedrals", "--bondangles"]
    if not len(fileTransmute) == len(list_of_valid_flags) and not len(fileTransmute) == (len(list_of_valid_flags) - 1):
        print_divide_between_command_and_output()
        print("Only five/six arguments are required, please check your input file.\n\n\n")

        options.print_help()
        sys.exit(0)

    for argument in fileTransmute:
        arg = argument.split()
        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
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

    return pdb_file, identifier, moiety, dihedrals, angles, conformation


def randomise(RandomiseInput, options):

    from transmute_func_tools import nucleoside_dict

    fileRandomise = remove_trailing_whitespace(list(map(lambda x: x.strip(), RandomiseInput.readlines())))

    # Check list of valid inputs 
    list_of_valid_flags = ["--chemistry", "--length", "--sequence", "--complement"]
    if len(fileRandomise) != 3:
        print_divide_between_command_and_output()
        print("Only three arguments are required at one time, please check our input file.\n\n\n")
        options.print_help()
        sys.exit(0)

    # Check if flags are valid
    for argument in fileRandomise:
        arg = argument.split()
        # If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
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

        if arg[0] == "--complement":
            if len(arg) > 2:
                complement = arg[1:]
                if check_if_nucleotides_are_valid(complement):
                    sys.exit(1)

            complement = arg[1]
            compl_test_against = ["homo"] + list(nucleoside_dict.keys())
            if complement not in compl_test_against:
                print("the variable in '--complement' is not a list and not one of the available inputs. Please revise your input for this flag.")
                sys.exit(0)

    # If the variable chemistry has not been defined, then
    try:
        chemistry
    except NameError:
        print_divide_between_command_and_output()
        print("No chemistry type was prompted! Revise your input file. \n")
        #options.print_help()
        sys.exit(0)

    return chemistry, length_sequence, sequence, complement


def xyz_to_pdb(ConversionInput, options):

    fileConversion = remove_trailing_whitespace(list(map(lambda x: x.strip(), ConversionInput.readlines())))

    # Check list of valid inputs 
    list_of_valid_flags = ["--xyz", "--atomID", "--atomname_list"]
    if len(fileConversion) != 3:
        print_divide_between_command_and_output()
        print("Only three arguments are required at one time. Please check your input file. \n\n\n ")
        options.print_help()
        sys.exit(0)

    # Check if flags are valid
    for argument in fileConversion:
        arg = argument.split()
        # If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
            #options.print_help()
            sys.exit(0)

        # Save the prompted arguments according to their respective flags
        if arg[0] == "--xyz":
            fname_xyz = arg[1]
            try:
                os.path.isfile(fname_xyz)
            except FileNotFoundError:
                print_divide_between_command_and_output()
                print(f"File {fname_xyz} was not found. Please revise either its name or its location on your system.\n")
                sys.exit(0)

        if arg[0] == "--atomID":
            atomID = arg[1]
            if len(atomID) > 4:
                print_divide_between_command_and_output()
                print("The argument '--atomID' requires up to three characters in its name.\n"
                        "Check your input of the atom identifier.\n")
                sys.exit(0)

        if arg[0] == "--atomname_list":
            atomname_list = arg[1:]

    return fname_xyz, atomID, atomname_list

