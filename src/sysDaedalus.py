from sys import version_info, exit
import os



""" This file exists to retrieve systems data on the current machine ."""

def version_checker():
    """ Checks the version of Daedalus"""
    if not version_info.major == 3 and version_info.minor >= 6:

        print("Python 3.6 or higher is required.")
        exit(1)

def return_DAEDALUS_home():
    """ Return the home directory of Daedalus """
    DAEDALUS_HOME = os.path.dirname(os.path.realpath(__file__)) + "/"

    return DAEDALUS_HOME

def exit_Daedalus():
    """ Exit Daedalus with the following message : """
    exit("\nExitted Daedalus unsuccesfully ... \n\n")


