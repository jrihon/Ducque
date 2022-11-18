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
    exit("\nExitted Ducque unsuccesfully ... \n\n")


