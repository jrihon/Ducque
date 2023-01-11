from sys import version_info, exit
import os



""" This file exists to retrieve systems data on the current machine ."""

def version_checker():
    """ Checks the version of Ducque"""
    if not version_info.major == 3 and version_info.minor >= 8:

        print("Python 3.8 or higher is required.")
        exit(1)

def return_DUCQUEHOME():
    """ Return the home directory of Ducque """
    DUCQUEHOME = os.path.dirname(
                        os.path.dirname(os.path.realpath(__file__))     #retrieve path of called module
                                ) + "/"

    return DUCQUEHOME

def exit_Ducque():
    """ Exit Ducque with the following message : """
    exit("\033[1;31m[ Exitted Ducque unsuccesfully ... ]\033[39m")



def print_inputfile(module, extensions):
    print(f"[INVALID INPUTFILE]: {module} input-files should end in {extensions}.")
    exit_Ducque()

def print_invalid_flag(flag : str):
    print(f"[INVALID FLAG]     : {flag}.")
    exit_Ducque()

def print_invalid_argument(argument : str, flag : str):
    print(f"[INVALID ARGUMENT] : {flag} does not receive `{argument}`.")
    exit_Ducque()

def print_insufficient_flag(amount : int):
    print(f"[INVALID QUERY]    : Insufficient amount of arguments. Required amount : {amount}.")
    exit_Ducque()


def check_for_empty_string(string : str) -> bool:

    try :
        string.split(sep=" ", maxsplit=1)[1]
    except IndexError : 
        print("Revise inputs!")
        exit_Ducque()

    else :
        s = string.split(sep=" ", maxsplit=1)[1]
        s = "".join(string) # join the entire list to one string

        if len(s.strip()) == 0 :
            return False
        return True


def print_build():                  print(f"[BUILD]           : Generating duplex.")
def print_time(t1, t0):             print(f"[TIME]            : %.4f seconds." % (t1 - t0))
def print_transmute(PDB_FNAME) :    print(f"[TRANSMUTE]       : Converting {PDB_FNAME} to a json file.")
def print_xyz(XYZ_FNAME):           print(f"[XYZ -> PDB]      : Converting {XYZ_FNAME} to a .pdb format file.")
def print_rand():                   print(f"[RANDOMISE]       : Randomisation of the given inputs!")
def print_dupl_ln(num_nucl):        print(f"[DUPLEX LENGTH]   : {num_nucl} nucleotides.")
def print_writing(outfile):         print(f"[WRITE FILE]      : {outfile} .")
def print_invalid_chemistry(chem):  print(f"[INVALID QUERY]   : `{chem}`. Please revise your inputs.")
def print_invalid_nucleoside(NA):   print(f"[INVALID QUERY]   : One or more of the nucleotides in the given sequence is invalid `{NA}`")
def print_empty_querty() :          print(f"[EMPTY QUERY]     : revisit your inputs")
def print_cant_find_Ducque():       print(f"[NOT FOUND]       : Ducque not found in the $PATH. Please add `Ducque` to the search path.")
