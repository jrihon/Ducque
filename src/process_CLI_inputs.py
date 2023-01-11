from builder import builder_library as LIB
import systemsDucque as SD
import os

""" When the user prompts the wrong values or flags, this python script intercepts most errors that happen at the start """

class InputExclusivity(Exception):
    Except = " These flags are mutually exclusive; --transmute     --build "

def remove_blank_lines(fileList : list) -> list:

    # intialise new list
    newList = []
    # Remove all trailing whitespace, cannot depend on only the last line being whitespace
    for i in fileList:
        if len(i.strip()) > 0:
            newList.append(i)

    return newList

def check_if_nucleotides_are_valid(input_sequence : list) -> bool:
    """ Check if any of the prompted nucleotides is not valid. """
    # Retrieve the keys of the dictionary from which we parse the data
    keysOfDict = LIB.codex_acidum_nucleicum.keys()
    # Check if any of the prompted nucleotides is not found in the sequence
    for NA in input_sequence:
        if NA not in keysOfDict:
            SD.print_invalid_nucleoside(NA)
            return False

    return True


def check_if_chemistry_is_valid(chemistry : str) -> str:
    """ If the chemistry is invalid, stop the script.
        If the chemistry is actually valid, return the keys of the dictionary as a list. """

    from transmute.transmute_constants import nucleoside_dict

    # Get the value from the transmute_constants nucleoside dictionary, to get the full name of the chemistry
    try:
        NUC_ID = nucleoside_dict[chemistry.upper()]
    except KeyError:
        SD.print_invalid_chemistry(chemistry)
        SD.exit_Ducque()

    return NUC_ID


def build(BUILDINPUT):

    if not BUILDINPUT.name.endswith(".binp"): 
        SD.print_inputfile("`--build`", "`.binp`")


    list_of_valid_flags = ["--sequence", "--complement", "--out"]
    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileDucque = remove_blank_lines(list(map(lambda x: x.strip(), BUILDINPUT.readlines())))

    # Check if amount of inputs are valid
    if len(fileDucque) != len(list_of_valid_flags):
        SD.print_insufficient_flag(3)

    for argument in fileDucque:
        flag = argument.split()[0]
        args = argument.split() # split the string into separate values

        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not flag in list_of_valid_flags:
            SD.print_invalid_flag(flag)

        # See if all arguments have been prompted
        SD.check_for_empty_string(argument)

        if flag == "--sequence":
            nucleicAcidList = list(map(lambda x: x.strip(","), args[1:]))

        if flag == "--complement":
            try:
                args[2] # if there is a third element, we are dealing with a list
            except:
                complement = args[1] # else it is just a single keyword
            else:
                complement =  list(map(lambda x: x.strip(","), args[1:]))

        if flag == "--out" :
            outFile = args[1]
    # If the variable outFile has not been prompted by the user, we default it ourselves by having it take on the name of the input file
    try :
        outFile
    except :
        outFile = BUILDINPUT.name.split(".")[0]
        list_of_valid_flags.remove("--out")

    # If the leading strand is valid
    if not check_if_nucleotides_are_valid(nucleicAcidList):
        SD.exit_Ducque()

    # if the complementary strand is valid
    if isinstance(complement, list) and len(complement) > 2 :
        if not check_if_nucleotides_are_valid(complement):
            SD.exit_Ducque()

    return nucleicAcidList, complement, outFile



def transmute(TRANSMUTEINPUT):

    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileTransmute = remove_blank_lines(list(map(lambda x: x.strip(), TRANSMUTEINPUT.readlines())))

    # Check if amount of inputs are valid
    list_of_valid_flags = ["--pdb", "--chemistry", "--moiety", "--conformation", "--dihedrals", "--bondangles"]
    if not len(fileTransmute) == len(list_of_valid_flags):
#    if not len(fileTransmute) == len(list_of_valid_flags) and not len(fileTransmute) == (len(list_of_valid_flags) - 1):
        SD.print_insufficient_flag(6)

    for argument in fileTransmute:
        flag = argument.split()[0]
        args = argument.split() 

        # Check if flags are valid. If a given flag is not a valid one, shut it down
        if not flag in list_of_valid_flags:
            SD.print_invalid_flag(flag)

        # See if all arguments have been prompted
        SD.check_for_empty_string(argument)

        # Save the prompted arguments according to the respective flags
        if flag == "--pdb":
            pdb_fname = args[1]

        if flag == "--chemistry":
            chemistry = args[1]

        if flag == "--conformation":
            conformation = args[1]

        if flag == "--moiety":
            moiety = args[1]

        if flag == "--dihedrals":
            dihedrals = args[1:]

        if flag == "--bondangles":
            angles = args[1:]

    return pdb_fname, chemistry, moiety, dihedrals, angles, conformation




def randomise(RANDOMISEINPUT):

    from transmute.transmute_constants import nucleoside_dict
    
    if not RANDOMISEINPUT.name.endswith(".rinp"): 
        SD.print_inputfile("`--randomise`", "`.rinp`")

    fileRandomise = remove_blank_lines(list(map(lambda x: x.strip(), RANDOMISEINPUT.readlines())))

    # Check list of valid inputs 
    list_of_valid_flags = ["--chemistry", "--length", "--sequence", "--complement"]
    if len(fileRandomise) != 3:
        SD.print_insufficient_flag(3)

    # Check if flags are valid
    for argument in fileRandomise:
        flag = argument.split()[0]
        args = argument.split()

        # If a given flag is not a valid one, shut it down
        if not flag in list_of_valid_flags:
            SD.print_invalid_flag(flag)

        # See if all arguments have been prompted
        SD.check_for_empty_string(argument)

        # Save the prompted arguments according to their respective flags
        if flag == "--chemistry":
            # If multiple chemistries are prompted
            if len(flag) > 2:
                chemistry = args[1:]
            # If one chemistry is prompted
            elif len(flag) == 2:
                chemistry = args[1]

        if flag == "--length":
            try : 
                int(args[1])
            except :
                print("Argument for `--length` is not an integer type. Try again.")
            else :
                length_sequence = int(args[1])
                sequence = None

        if flag == "--sequence":
            sequence = args[1:]
            length_sequence = 0

        if flag == "--complement":
            if len(args) > 2:
                complement =  list(map(lambda x: x.strip(","), args[1:]))
                if not check_if_nucleotides_are_valid(complement):
                    SD.exit_Ducque()
                continue

            # Check if the complement flag contains either the 'homo' flag or a valid chemistry.
            complement = args[1]
            compl_test_against = ["homo"] + list(nucleoside_dict.keys())
            if complement not in compl_test_against:
                SD.print_invalid_argument(complement, "--complement")

#    # If the variable chemistry has not been defined, then
#    try:
#        chemistry
#    except NameError:
#        SD.print_invalid_argument('', "--chemistry")

    outFile = RANDOMISEINPUT.name.split(".")[0]
    return chemistry, length_sequence, sequence, complement, outFile


def xyz_to_pdb(CONVERSIONINPUT):

    fileConversion = remove_blank_lines(list(map(lambda x: x.strip(), CONVERSIONINPUT.readlines())))

    # Check list of valid inputs 
    list_of_valid_flags = ["--xyz", "--residue", "--atomname_list"]
    if len(fileConversion) != 3:
        SD.print_insufficient_flag(3)

    # Check if flags are valid
    for argument in fileConversion:
        flag = argument.split()[0]
        args = argument.split()
        # If a given flag is not a valid one, shut it down
        if not flag in list_of_valid_flags:
            SD.print_invalid_flag(flag)

        # See if all arguments have been prompted
        SD.check_for_empty_string(argument)

        # Save the prompted arguments according to their respective flags
        if flag == "--xyz":
            xyz_fname = args[1]
            try:
                os.path.isfile(xyz_fname)
            except FileNotFoundError:
                print(f"File {xyz_fname} was not found. Please revise either its name or its location on your system.\n")
                SD.exit_Ducque()

        if flag == "--residue":
            atomID = args[1]
            if len(atomID) > 4:
                print("The argument '--atomID' requires up to three characters in its name.\n"
                        "Check your input of the atom identifier.\n")
                SD.exit_Ducque()

        if flag == "--atomname_list":
            atomname_list = args[1:]

    return xyz_fname, atomID, atomname_list

def gui_module(GUIINPUT): 

    # Take GUIINPUT out of the list
    GUIINPUT = GUIINPUT[0]

    if GUIINPUT == "NO_FLAG" : 
        return GUIINPUT

    list_of_valid_flags = ["build", "transmute", "randomise", "xyz_pdb"]

    if not GUIINPUT in list_of_valid_flags:
        SD.print_invalid_flag(GUIINPUT)

    return GUIINPUT
