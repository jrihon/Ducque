import labyrinth_repository as LabRepo
import sysDaedalus
import os

""" When the user prompts the wrong values or flags, this python script intercepts most errors that happen at the start """

class InputExclusivity(Exception):
    Except = " These flags are mutually exclusive; --transmute     --Daedalus "




def print_divide_between_command_and_output():
    """ This function exists solely to split the Daedalus call command and its output."""
    print("-----------------------------------------------------------")


def remove_blank_lines(fileList : list) -> list:

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
    keysOfDict = LabRepo.codex_acidum_nucleicum.keys()
    # Check if any of the prompted nucleotides is not found in the sequence

    for NA in input_sequence:
        if NA not in keysOfDict:
            print_divide_between_command_and_output()
            print(f"One or more of the nucleotides in the given sequence is invalid. Please check your input file : {NA} \n\n")
            return False

    return True


def check_if_chemistry_is_valid(chemistry : str) -> str:
    """ If the chemistry is invalid, stop the script.
        If the chemistry is actually valid, return the keys of the dictionary as a list. """

    from transmute_func_tools import nucleoside_dict

    # Get the value from the transmute_func_tools nucleoside dictionary, to get the full name of the chemistry
    try:
        NUC_ID = nucleoside_dict[chemistry.upper()]
    except KeyError:
        print(f"The following key does not exist in the dictionary : {chemistry.upper()}.\nPlease revise your inputs.\n")
        sysDaedalus.exit_Daedalus()

    return NUC_ID


def daedalus(DAEDALUSINPUT, options):
    print("-----------------------------------------------------------")
    print("Daedalus - Nucleic Acid Architecture initiated! Building sequence ...\n")

    list_of_valid_flags = ["--sequence", "--complement", "--out"]
    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileDaedalus = remove_blank_lines(list(map(lambda x: x.strip(), DAEDALUSINPUT.readlines())))

    for argument in fileDaedalus:
        arg = argument.split()

        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
            options.print_help()
#            sysDaedalus.exit_Daedalus()

        if arg[0] == "--sequence":
            nucleicAcidList = list(map(lambda x: x.strip(","), arg[1:]))

        if arg[0] == "--complement":
            # If there is a input possibility at index 2, meaning more than one string have been inputted, then get the entire string as a list variable.
            # If there is not an input possibility at index 2, this means there is only one input after the flag available and that means it is just a string.
            try:
                arg[1]
            except:
                print_divide_between_command_and_output()
                print("You have forgotten to include an argument for the --complement flag! Please add this to your input file.")
                sysDaedalus.exit_Daedalus()

            try:
                arg[2]
            except:
                complement = arg[1]
            else:
                complement =  list(map(lambda x: x.strip(","), arg[1:]))

        if arg[0] == "--out" :
            outFile = arg[1]

    # If the variable outFile has not been prompted by the user, we default it ourselves by having it take on the name of the input file
    try :
        outFile
    except :
        outFile = DAEDALUSINPUT.name.split(".")[0]
        list_of_valid_flags.remove("--out")

    # Check if amount of inputs are valid
    if len(fileDaedalus) != len(list_of_valid_flags):
        print_divide_between_command_and_output()
        print("Only two/three arguments are required, please check our input file.\n`--out` is an optional flag.\n\n\n")
#        options.print_help()
        sysDaedalus.exit_Daedalus()

    # If a given nucleotide is not in the list of valid nucleotides, stop the program
    if not check_if_nucleotides_are_valid(nucleicAcidList):
        #options.print_help()
        sysDaedalus.exit_Daedalus()

    if isinstance(complement, list) and len(complement) > 2 :
        if not check_if_nucleotides_are_valid(complement):
            #options.print_help()
            sysDaedalus.exit_Daedalus()

    return nucleicAcidList, complement, outFile


def transmute(TRANSMUTEINPUT, options):

    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileTransmute = remove_blank_lines(list(map(lambda x: x.strip(), TRANSMUTEINPUT.readlines())))

    # Check if amount of inputs are valid
    list_of_valid_flags = ["--pdb", "--chemistry", "--moiety", "--conformation", "--dihedrals", "--bondangles"]
    if not len(fileTransmute) == len(list_of_valid_flags) and not len(fileTransmute) == (len(list_of_valid_flags) - 1):
        print_divide_between_command_and_output()
        print("Only five/six arguments are required, please check your input file.\n`--conformation` is an optional flag.\n\n\n")
        options.print_help()
#        sysDaedalus.exit_Daedalus()

    for argument in fileTransmute:
        arg = argument.split()
        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
            #options.print_help()
            sysDaedalus.exit_Daedalus()

        # Save the prompted arguments according to the respective flags
        if arg[0] == "--pdb":
            pdb_fname = arg[1]

        if arg[0] == "--chemistry":
            chemistry = arg[1]

        if arg[0] == "--conformation":
            conformation = arg[1]

        if arg[0] == "--moiety":
            moiety = arg[1]

        if arg[0] == "--dihedrals":
            dihedrals = arg[1:]

        if arg[0] == "--bondangles":
            angles = arg[1:]

        # if the conformation was not prompted, since it's an optional argument
        try:
            conformation
        except NameError:
            conformation = False

    return pdb_fname, chemistry, moiety, dihedrals, angles, conformation


def randomise(RANDOMISEINPUT, options):

    from transmute_func_tools import nucleoside_dict

    fileRandomise = remove_blank_lines(list(map(lambda x: x.strip(), RANDOMISEINPUT.readlines())))

    # Check list of valid inputs 
    list_of_valid_flags = ["--chemistry", "--length", "--sequence", "--complement"]
    if len(fileRandomise) != 3:
        print_divide_between_command_and_output()
        print("Only three arguments are required at one time, please check our input file.\n\n\n")
        options.print_help()
#        sysDaedalus.exit_Daedalus()

    # Check if flags are valid
    for argument in fileRandomise:
        arg = argument.split()
        # If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
            #options.print_help()
            sysDaedalus.exit_Daedalus()

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
                complement =  list(map(lambda x: x.strip(","), arg[1:]))
                if not check_if_nucleotides_are_valid(complement):
                    sysDaedalus.exit_Daedalus()
                continue

            # Check if the complement flag contains either the 'homo' flag or a valid chemistry.
            complement = arg[1]
            compl_test_against = ["homo"] + list(nucleoside_dict.keys())
            if complement not in compl_test_against:
                print("the input of '--complement' is not a valid list and not one of the possible inputs. Please revise your input for this flag.")
                sysDaedalus.exit_Daedalus()

    # If the variable chemistry has not been defined, then
    try:
        chemistry
    except NameError:
        print_divide_between_command_and_output()
        print("No chemistry type was prompted! Revise your input file. \n")
        #options.print_help()
        sysDaedalus.exit_Daedalus()

    outFile = RANDOMISEINPUT.name.split(".")[0]
    return chemistry, length_sequence, sequence, complement, outFile


def xyz_to_pdb(CONVERSIONINPUT, options):

    fileConversion = remove_blank_lines(list(map(lambda x: x.strip(), CONVERSIONINPUT.readlines())))

    # Check list of valid inputs 
    list_of_valid_flags = ["--xyz", "--atomID", "--atomname_list"]
    if len(fileConversion) != 3:
        print_divide_between_command_and_output()
        print("Only three arguments are required at one time. Please check your input file. \n\n\n ")
        options.print_help()
#        sysDaedalus.exit_Daedalus()

    # Check if flags are valid
    for argument in fileConversion:
        arg = argument.split()
        # If a given flag is not a valid one, shut it down
        if not arg[0] in list_of_valid_flags:
            print_divide_between_command_and_output()
            print(f"\n\nThe following flag is invalid : {arg[0]}. Please check your input file.\n\n\n")
            #options.print_help()
            sysDaedalus.exit_Daedalus()

        # Save the prompted arguments according to their respective flags
        if arg[0] == "--xyz":
            xyz_fname = arg[1]
            try:
                os.path.isfile(xyz_fname)
            except FileNotFoundError:
                print_divide_between_command_and_output()
                print(f"File {xyz_fname} was not found. Please revise either its name or its location on your system.\n")
                sysDaedalus.exit_Daedalus()

        if arg[0] == "--atomID":
            atomID = arg[1]
            if len(atomID) > 4:
                print_divide_between_command_and_output()
                print("The argument '--atomID' requires up to three characters in its name.\n"
                        "Check your input of the atom identifier.\n")
                sysDaedalus.exit_Daedalus()

        if arg[0] == "--atomname_list":
            atomname_list = arg[1:]

    return xyz_fname, atomID, atomname_list

