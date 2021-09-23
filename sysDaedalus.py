import sys
import os



""" This file exists to retrieve systems data on the current machine ."""

def version_checker():
    """ Checks the version of Daedalus"""
    import sys
    if not sys.version_info.major == 3 and sys.version_info.minor >= 6:

        print("Python 3.6 or higher is required.")
        sys.exit(1)

def return_DAEDALUS_home():
    """ Return the home directory of Daedalus """
    daedalus_home = os.path.dirname(os.path.realpath(__file__)) + "/"

    return daedalus_home

