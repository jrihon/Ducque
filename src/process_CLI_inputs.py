import systemsDucque as SD
import os

from transmute.utils_transmute import LinkerTransmute, NucleosideTransmute, NucleobaseTransmute
from ducquelib.library import TABLE_NUCLEOTIDES, TABLE_CHEMISTRY
from nucleobase.utils_nucleobase import NucleobaseMod

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
    keysOfDict = TABLE_NUCLEOTIDES.keys()
    # Check if any of the prompted nucleotides is not found in the sequence
    for NA in input_sequence:
        if NA not in keysOfDict:
            SD.print_invalid_nucleoside(NA)
            return False

    return True


def check_if_chemistry_is_valid(chemistry : str) -> str:
    """ If the chemistry is invalid, stop the script.
        If the chemistry is actually valid, return the keys of the dictionary as a list. """


    # Get the value from the transmute_constants nucleoside dictionary, to get the full name of the chemistry
    try:
        NUC_ID = TABLE_CHEMISTRY[chemistry.upper()]
    except KeyError:
        SD.print_invalid_chemistry(chemistry)
        SD.exit_Ducque()

    return NUC_ID


def build(BUILDINPUT):

    if not BUILDINPUT.name.endswith(".binp"): 
        SD.print_inputfile("`--build`", "`.binp`")


    list_of_valid_flags = ["--sequence", "--complement", "--pdbname"]
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
        SD.check_for_empty_string(argument, flag)

        if flag == "--sequence":
            nucleicAcidList = list(map(lambda x: x.strip(",").upper(), args[1:]))

        if flag == "--complement":
            try:
                args[2] # if there is a third element, we are dealing with a list
            except:
                complement = args[1] # else it is just a single keyword
            else:
                complement =  list(map(lambda x: x.strip(",").upper(), args[1:]))

        if flag == "--pdbname" :
            pdbName = args[1]


    # If the leading strand is valid
    if not check_if_nucleotides_are_valid(nucleicAcidList):
        SD.exit_Ducque()

    # if the complementary strand is valid
    if isinstance(complement, list) and len(complement) > 2 :
        if not check_if_nucleotides_are_valid(complement):
            SD.exit_Ducque()

    return nucleicAcidList, complement, pdbName



def transmute(TRANSMUTEINPUT) -> LinkerTransmute | NucleosideTransmute | NucleobaseTransmute :

    if not TRANSMUTEINPUT.name.endswith(".tinp"): 
        SD.print_inputfile("`--transmute`", "`.tinp`")

    # Read the input file and create a list object of the sequence. Removes any whitespace.
    fileTransmute = remove_blank_lines(list(map(lambda x: x.strip(), TRANSMUTEINPUT.readlines())))

#    # Check if amount of inputs are valid
#    list_of_valid_flags = ["--pdb", "--chemistry", "--moiety", "--conformation", "--dihedrals", "--bondangles"]
    # iterate over the possible inputs to check if `nucleoside`, `linker` or `nucleobase` is prompted
    moiety = ""
    for x in fileTransmute:
        line = x.split() 
        flag = line[0]
        if flag == "--moiety": 
            SD.check_for_empty_string(x, flag)
            if line[1] == "nucleoside": moiety += "nucleoside"
            elif line[1] == "linker": moiety += "linker"
            elif line[1] == "nucleobase": moiety += "nucleobase"
            else : 
                SD.print_invalid_argument("--moiety", line[1])

    # Prepare nucleoside transmutation
    if moiety == "nucleoside" : 

        nucT = NucleosideTransmute()
        nucT.set_attributes(fileTransmute)

        return nucT

    # Prepare linker transmutation
    elif moiety == "linker" : 

        linkT = LinkerTransmute()
        linkT.set_attributes(fileTransmute)

        return linkT


    # Prepare nucleobase transmutation
    elif moiety == "nucleobase" : 

        baseT = NucleobaseTransmute()
        baseT.set_attributes(fileTransmute)

        return baseT


    else : 
        SD.print_empty_query("--moiety")

#    if not len(fileTransmute) == len(list_of_valid_flags):
#        SD.print_insufficient_flag(len(list_of_valid_flags))
#
#    for argument in fileTransmute:
#        flag = argument.split()[0]
#        args = argument.split() 
#
#        # Check if flags are valid. If a given flag is not a valid one, shut it down
#        if not flag in list_of_valid_flags:
#            SD.print_invalid_flag(flag)
#
#        # See if all arguments have been prompted
#        SD.check_for_empty_string(argument, flag)
#
#        # Save the prompted arguments according to the respective flags
#        if flag == "--pdb":
#            pdb_fname = args[1]
#
#        if flag == "--chemistry":
#            chemistry = args[1]
#
#        if flag == "--conformation":
#            conformation = args[1]
#
#        if flag == "--moiety": 
#            moiety = args[1]
#
#        if flag == "--nucleobase": # `nucleoside only` flag
#            nucleobase = args[1]
#
#        if flag == "--dihedrals":
#            dihedrals = args[1:]
#
#        if flag == "--bondangles":
#            angles = args[1:]
#
#    # If nucleobase is not prompted, the linker moiety is queried, make this variable empty
#    try : 
#        nucleobase
#    except :
#        nucleobase = ""
#
#    return pdb_fname, chemistry, moiety, dihedrals, angles, conformation, nucleobase




def randomise(RANDOMISEINPUT):

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
        SD.check_for_empty_string(argument, flag)

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
                SD.print_invalid_argument(args[1], "--length")
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
            compl_test_against = ["HOMO"] + list(TABLE_CHEMISTRY.keys())
            if complement not in compl_test_against:
                SD.print_invalid_argument(complement, "--complement")

#    # If the variable chemistry has not been instantiated, then
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
        SD.check_for_empty_string(argument, flag)

        # Save the prompted arguments according to their respective flags
        if flag == "--xyz":
            xyz_fname = args[1]
            try:
                os.path.isfile(xyz_fname)
            except FileNotFoundError:
                SD.print_filenotfound(xyz_fname)
                SD.exit_Ducque()

        if flag == "--residue":
            residue = args[1]
            if len(residue) > 4:
                SD.print_invalid_argument(args[1], "--residue")

        if flag == "--atomname_list":
            atomname_list = args[1:]

    return xyz_fname, residue, atomname_list

def gui_module(GUIINPUT): 

    # Take GUIINPUT out of the list
    GUIINPUT = GUIINPUT[0]

    if GUIINPUT == "NO_FLAG" : 
        return GUIINPUT

    list_of_valid_flags = ["build", "transmute", "randomise", "xyz_pdb", "tlinker"]

    if not GUIINPUT in list_of_valid_flags:
        SD.print_invalid_flag(GUIINPUT)

    return GUIINPUT


def nucleobase(NBASEINPUT) -> tuple[str, list[NucleobaseMod]]:

    if not NBASEINPUT.name.endswith(".ninp"): 
        SD.print_inputfile("`--nbase`", "`.ninp`")

    fileNbases = remove_blank_lines(list(map(lambda x: x.strip(), NBASEINPUT.readlines())))
    list_of_modifications = list()

    pdb_nbase_fname = ""

    # find the --pdb flag 
    for line in fileNbases : 
        lList = line.split()
        if lList[0] == "--pdb" : 
            try : 
                pdb_nbase_fname = lList[1]
            except IndexError :
                SD.print_empty_query("--pdb")
                SD.exit_Ducque()

            if not os.path.isfile(pdb_nbase_fname): 
                SD.print_filenotfound(pdb_nbase_fname)
                SD.exit_Ducque()

            break

    # If the --pdb flag is not found, exit early
    if pdb_nbase_fname == '': 
        SD.print_empty_query("--pdb")
        SD.exit_Ducque()


    argumentless_flags = ["--nucleobase" ,"-reorient"]
    # Fill out `list_of_modifications` with NucleobaseMod objects
    isFirstObjectedCreated = False
    nb_instance=  0
    for argument in fileNbases : 
        flag = argument.split()[0]

        # nucleobase and reorient are two flags that do not take an argument
        # if we cannot index into a second element of the other flags, crash it
        if flag not in argumentless_flags: 
            try : 
                arg = argument.split()[1]
            except IndexError : 
                SD.print_empty_query(flag)
                SD.exit_Ducque()
            

        # if line starts with `--nucleobase`
        if flag == "--nucleobase" : 

            # if there are no entries yet, instance a NucleobaseMod()
            if not isFirstObjectedCreated : 
                nb_instance += 1
                nbase = NucleobaseMod(nb_instance)
                isFirstObjectedCreated = True
                continue

            else : 
                # If we encounter another instance of --nucleobase
                # append the current one to the list
                list_of_modifications.append(nbase)

                # Clean wipe data and initalise new NucleobaseMod() class
                nb_instance += 1
                nbase = NucleobaseMod(nb_instance)
                continue


        if flag == "-position": 
            nbase.set_position(arg, pdb_nbase_fname)

        if flag == "-mod": 
            nbase.set_mod(arg)

        if flag == "-resname": 
            nbase.set_new_resname(arg)

        if flag == "-reorient": 
            nbase.set_reorient()

    # try to append last instance of --nucleobase
    try : 
        list_of_modifications.append(nbase)
    # if it fails, that means the --nucleobase flag was never encountered and a NucleobaseMod() object was never instanced
    except : 
        SD.print_empty_query("--nucleobase")
        SD.exit_Ducque()
            
    # Final check before releasing it into the wild
    for nb in list_of_modifications : 
        nb.check_empty_attributes()                     # checks for empty position and empty mod arguments
        nb.check_if_resname_queried(pdb_nbase_fname)    # checks if residue name is passed in argument or 
                                                        #   can be parsed from given pdb

    return pdb_nbase_fname, list_of_modifications
