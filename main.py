import argparse
import sys
import os
from time import time

import transmute    # Convert pdb to json or vice versa
import labyrinth    # Build the duplex
import fundaments   # Exceptions or custom errors

explanation = """
Project to generate and customize DNA, RNA, XNA duplex molecules

Designed and written by Doctorandus Rihon Jérôme.\n
______________________________________________________________________
"""
t1 = time()
# ---------------------------------------- P A R S E  A R G U M E N T S ----------------------------------------
options = argparse.ArgumentParser(description=explanation,
                                  add_help=False,
                                  formatter_class=argparse.RawTextHelpFormatter)

## Let's not make this argument a required one, since we also might want to convert json to pdb
options.add_argument('--pdb',
                     help='Input the pdb of the nucleic acid you want to convert to json')

options.add_argument('--json',
                     help='Input the json of the nucleic you want to convert to pdb')

options.add_argument('--Daedalus',
        help='Ask the software what you want it to build')

options.add_argument('--help', action='help',
        help='Input the duplex of interest with the " -- Daedalus " flag.')

arguments = options.parse_args()

if len(sys.argv) == 1:
    options.print_help()
    sys.exit(0)

try:
    if arguments.pdb and arguments.json:
        raise fundaments.InputExclusivity

except fundaments.InputExclusivity:
    print("InputExclusivity: These flags are mutually exclusive; --pdb --json ")
    sys.exit(0)

# ---------------------------------------- M A I N ----------------------------------------
def main():

    # Convert pdb to json
    if arguments.pdb and transmute.pdb_is_linker(arguments.pdb):
        print('Reading linker molecule ...')
        transmute.linker_to_json(arguments.pdb)
        sys.exit(0)
    elif arguments.pdb:
        print('Reading nucleic acid molecule ...')
        transmute.NA_to_json(arguments.pdb)
        sys.exit(0)

    # Convert json to pdb
    if arguments.json and transmute.json_is_linker(arguments.json):
        print('Reading linker .json ...')
        transmute.linker_to_pdb(arguments.json)
        sys.exit()
    elif arguments.json:
        print('Reading nucleic acid .json ...')
        transmute.NA_to_pdb(arguments.json)
        sys.exit(0)

    if arguments.Daedalus:
        labyrinth.Architecture(arguments.Daedalus)

    print('Time spent: %.5f seconds.' % (time() - t1))
# ---------------------------------------- S T A R T   T H E   S H O W ----------------------------------------
if __name__ == '__main__':
    main()
