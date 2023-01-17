from numpy import asarray
import os
import sys


"""
    Generate a standard prep file with the input of a given pdb file and a file that contains the molecule tree structure and accompanying atomic charges

    First argument :
        a file that contains the following and this exact order --

        ´´´´´´´´´
        atom names      amber tree structure        atom typing         RESP charges
        atom names      amber tree structure        atom typing         RESP charges
        atom names      amber tree structure        atom typing         RESP charges
        atom names      amber tree structure        atom typing         RESP charges
        atom names      amber tree structure        atom typing         RESP charges
         ....               .....                       ....                ....
        atom names      amber tree structure        atom typing         RESP charges
        atom names      amber tree structure        atom typing         RESP charges
        ´´´´´´´´´


        The spacing in between does not matter, but there has to be whitespace in between the columns
        The quotation marks are there to just denote the contents on the file and are not part of the actual inputs

    Second argument :
        The name of the residue
        This should be either two (2) or three (3) characters long at best.
            (This is pdb convention, not my rule)


    The rest of the prep file should be filled in manually by the user themselves.
    Mind the IMPROPER and LOOP CLOSING EXPLICIT section, especially important in nucleic acids

    CAVEAT : Make sure there are not trailing whitespaces or trailing lines underneath the columns

    For more information, consult the `How_To_Create_A_Forcefield.pdf` manual
"""





def main():

    ## First Argument
    try:
        molecule_inputs = sys.argv[1]
    except FileNotFoundError:
        raise
    else:
        # read in data
        atomnames = list()
        tree = list()
        types = list()
        charges = list()

        with open(molecule_inputs, "r") as MI :
            for line in MI:
                data_tuple = line.split()
                atomnames.append(data_tuple[0])
                tree.append(data_tuple[1])
                types.append(data_tuple[2])
                charges.append(data_tuple[3])


    ## Second argument - residue
    try:
        residue_name = sys.argv[2]
    except IndexError:
        raise

    # strip any whitespace and empty strings from the list
    for i in atomnames, tree, types, charges :
        i = [x for x in list(map(lambda x : x.strip(), i)) if x]

    # Create a range of indexes and also the lists of 
    range_index = range_indexes(atomnames)

    bonds = range_index - 1
    angles = range_index - 2
    dihedrals = range_index - 3

    current_dir = os.path.basename(os.getcwd())
    with open(current_dir + ".prep", "w") as prep :


        # Write the header of the prep file
        prep.write(header_prep(residue_name))

        # Write out all the available data
        for i in range(len(range_index)):
            output_tuple = [range_index[i],
                            atomnames[i],
                            types[i],
                            tree[i],
                            bonds[i],
                            angles[i],
                            dihedrals[i],
                            1.20,
                            120.00,
                            180.00,
                            charges[i]]

            prep.write("%4d   %-4s  %-2s    %s %4d%4d%4d    %1.2f    %3.2f    %4.2f%13s\n" % tuple(output_tuple))


        # Write the last parts of the prep file
        prep.write(tail_prep())


def range_indexes(genericList : list) -> list:

    return asarray([x + 4 for x in range(len(genericList))])


def header_prep(residue : str) -> str:

    if len(residue) > 3 or len(residue) < 2:
        sys.exit("The amount of character in the residue should either be two or three characters long.")

    header ="SAMPLE TEXT - INFORMATION ON CHEMISTRY AND FRAGMENT\n\n" \
    " " + residue + "  INT     1\n" \
    " CORRECT NOMIT DU   BEG\n" \
    "   0.000\n" \
    "   1   DUMM  DU    M    0  -1  -2    0.00      0.00      0.00       0.0000\n" \
    "   2   DUMM  DU    M    1   0  -1    1.00      0.00      0.00       0.0000\n" \
    "   3   DUMM  DU    M    2   1   0    1.00     90.00      0.00       0.0000\n"


    return header

def tail_prep() -> str:

    tail ="\n" \
    "IMPROPER\n" \
    "\n" \
    "LOOP CLOSING EXPLICIT\n" \
    "\n" \
    "DONE\n"

    return tail



main()
