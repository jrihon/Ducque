"""
The python file that organises the labyrinth_func.py better
"""
import numpy as np
import json

import labyrinth

def check_slope_of_array(arr : np.array) -> str:
    """ In labyrinth_func.py there is a function that retrieves the interpolated dihedral angle
    But it works on whether or not the list is ascending or descending
    That's what we need to figure out here now and return this """

    slope_list = []
    for i in range(len(arr)):
        if i == (len(arr) - 1):
            if arr[i] > arr[0]:
                slope_list.append("D")
            else:
                slope_list.append("A")

        elif arr[i] > arr[i+1]:
            slope_list.append("D")
        else:
            slope_list.append("A")

    descending = slope_list.count("D")
    ascending = slope_list.count("A")

    if ascending > descending:
        return "ASCENDING"
    else:
        return "DESCENDING"


def retrieve_atom_index(json_object, atom : str) -> np.array :
    """ Retrieves the index in the json_object.array of the atom of interest
    This integer will be used to retrieve the vector of the atom of interest """
    return json_object.atom_list.index(atom)


def pdb_AtomNames_or_ElementSymbol(list_of_sequence : list, identifier : str) -> list:
    """ Loads in the atom names from the json files
    The identifier is either the string "Atoms" or the string "ElementSymbol", which will parse the list of interest """

    # reverse the sequence of the list, because we build to sequence from the bottom up
    list_of_sequence.reverse()

    # Initialise an empty list
    atom_list = []

    for i in range(len(list_of_sequence)):

        # If it is the first in the nucleotide sequence, meaning the last in the nucleotide sequence
        if i == 0:
            buildingblock = list_of_sequence[i]
            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            tmp_atomlist = json.loads(nucleoside["pdb_properties"][identifier])
            atom_list = atom_list + tmp_atomlist

        # If it is the last in the reversed sequence, meaning the first nucleotide
        if (i + 1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = labyrinth.codex_acidum_nucleicum[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            tmp_atomlist = json.loads(linker["pdb_properties"][identifier]) + json.loads(nucleoside["pdb_properties"][identifier])

            atom_list = atom_list + tmp_atomlist

            return atom_list

        if (i + 1) == len(list_of_sequence):
            # since this is not the last one, just carry on as usual
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = labyrinth.codex_acidum_nucleicum[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            tmp_atomlist = json.loads(linker["pdb_properties"][identifier]) + json.loads(nucleoside["pdb_properties"][identifier])

            atom_list = atom_list + tmp_atomlist


def pdb_Sequence(list_of_sequence : list) -> np.array :
    """ Determines the number in the sequence of the nucleotides in the strands based off on the shape of their array  """

    # reverse the sequence of the list, because we build to sequence from the bottom up
    list_of_sequence.reverse()

    # Initialise an empty array
    sequence_array = np.array([], dtype=int)

    # Initialise a counter
    seq_count = 1

    for i in range(len(list_of_sequence)):

        # If it is the first in the nucleotide sequence, meaning the last in the nucleotide sequence
        if i == 0:
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            nuc_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0]
            tmp_seqarray = np.full(nuc_shape, seq_count, dtype=int)

            seq_count += 1

            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)

        # If it is the last in the reversed sequence, meaning the first nucleotide
        if (i + 1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = labyrinth.codex_acidum_nucleicum[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
            tmp_seqarray = np.full(nucleotide_shape, seq_count, dtype=int)

            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)

            # Return the output of the atom_list, as this is the last nucleotide in the sequence
            return sequence_array

        # since this is not the last one or the first one, just carry on as usual
        if not i == 0 and not (i+1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = labyrinth.codex_acidum_nucleicum[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
            tmp_seqarray = np.full(nucleotide_shape, seq_count, dtype=int)

            seq_count += 1

            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)


def pdb_Residuename(list_of_sequence : list) -> list:

    # reverse the sequence of the list, because we build to sequence from the bottom up
    list_of_sequence.reverse()

    # Initialise an empty array
    resname_list = []

    for i in range(len(list_of_sequence)):

        # If it is the first in the nucleotide sequence, meaning the last in the nucleotide sequence
        if i == 0:
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            nuc_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0]
            ID = json.loads(nucleoside["Identity"])[2]

            tmp_resname = [ID for i in range(nuc_shape)]
            resname_list = resname_list + tmp_resname

        # If it is the last nucleotide in the sequence
        if (i + 1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = labyrinth.codex_acidum_nucleicum[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
            ID = json.loads(nucleoside["Identity"])[2]

            tmp_resname = [ID for i in range(nucleotide_shape)]
            resname_list = resname_list + tmp_resname
            # Return the output of the atom_list, as this is the last nucleotide in the sequence
            return resname_list


        # since this is not the last one or the first one, just carry on as usual
        if not i == 0 and not (i+1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = labyrinth.codex_acidum_nucleicum[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = labyrinth.codex_acidum_nucleicum[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
            ID = json.loads(nucleoside["Identity"])[2]

            tmp_resname = [ID for i in range(nucleotide_shape)]
            resname_list = resname_list + tmp_resname

