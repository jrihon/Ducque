from sys import version_info, exit
import os



""" 
    This file exists to retrieve systems data on the current machine .
    
    This also exists to return error messages from I/O operations.
"""

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
    # prints the text in red colors and follows it up with a revert to the original terminal text color (white)



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

def angle_exclusivity():
    print("[EXCL. ANGLES]    : Cannot have η-dihedral without ζ-dihedral!")
    print("[REFERENCE]       : 1983 IUPAC on Nucleic Acids : ` 1982 Mar 1;131(1):9-15. doi: 10.1111/j.1432-1033.1983.tb07225.x.`")
    print("[LINK]            : https://pubmed.ncbi.nlm.nih.gov/6832147/")

def print_no_overwrite(fname, d):
    print(f"\033[34m[OVERWRITE BLOCK] : Cannot overwrite. File already present in Current Directory ")
    print(f"    [FILE]        : `{fname}`")
    print(f"    [DIRECTORY]   : `{d}`\033[39m")

def print_time(t1, t0):             print(f"[TIME]            : %.4f seconds." % (t1 - t0))
def print_build():                  print(f"\033[96m[BUILD]           : Generating duplex.\033[39m\n")
def print_transmute(PDB_FNAME) :    print(f"\033[96m[TRANSMUTE]       : Converting {PDB_FNAME} to a json file.\033[39m\n")
def print_xyz(XYZ_FNAME):           print(f"\033[96m[XYZ -> PDB]      : Converting {XYZ_FNAME} to a .pdb format file.\033[39m\n")
def print_rand():                   print(f"\033[96m[RANDOMISE]       : Randomisation of the given inputs!\033[39m\n")
def print_dupl_ln(num_nucl):        print(f"[DUPLEX LENGTH]   : {num_nucl} nucleotides.")
def print_writing(outfile):         print(f"\033[94m[WRITE FILE]      : {outfile} .\033[39m")
def print_invalid_chemistry(chem):  print(f"[INVALID QUERY]   : `{chem}`. Please revise your inputs.")
def print_invalid_nucleoside(NA):   print(f"[INVALID QUERY]   : One or more of the nucleotides in the given sequence is invalid `{NA}`")
def print_empty_query(flag) :       print(f"[EMPTY QUERY]     : No input found for `{flag}`.")
def print_cant_find_Ducque():       print(f"[NOT FOUND]       : Ducque not found in the $PATH. Please add `Ducque` to the search path.")
def print_insuf_amount(flag):       print(f"[INVALID AMOUNT]  : Incorrect amount of queries to properly fill the `{flag}` entry")
def print_conversion_err(name, x):  print(f"[CONVERSION ERROR]: Could not convert the angle `{name}` to a float `{x}` ")
def print_launch(module):           print(f"[LAUNCH MODULE]   : {module.title()} module !")
def print_filenotfound(fname):      print(f"[FILE NOT FOUND]  : The following file was not found at the current path `{fname}`.")
