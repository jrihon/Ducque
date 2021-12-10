import numpy as np
import json
import os
from typing import Union

import labyrinth
import labyrinth_func_tools3 as LFT3    # Parse the nucleic acid dictionaries
""" This scripts makes data parsing much easier and makes labyrinth.py much more organised. """

CODEX = LFT3.codex_acidum_nucleicum

def check_slope_of_array(arr : np.array) -> str:
    """ In labyrinth_func_tools1.py there is a function that retrieves the interpolated dihedral angle
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


def retrieve_atom_index(json_object, atom : str, index_counter : int = 0) -> int :
    """ Retrieves the index in the json_object.array of the atom of interest
    This integer will be used to retrieve the vector of the atom of interest """
    return json_object.atom_list.index(atom) + index_counter


def retrieve_atom_index_MULTIPLE(json_object, atoms : list, index_counter : int = 0) -> np.array :
    """ Retrieves the index in the json_object.array of the atom of interest
        This integer will be used to retrieve the vector of the atom of interest """
    array_of_indexes = np.zeros(len(atoms), dtype=int)

    for i in range(len(atoms)):
        array_of_indexes[i] = json_object.atom_list.index(atoms[i]) + index_counter

    return array_of_indexes


def retrieve_list_of_dihedrals_and_angles_to_build_with(next_nucleoside) -> list:
    """ Retrieve a list of the names of the angles and dihedrals we need to parse from the jsonObjects to further build our nucleic acid strand."""
    # Retrieve all the dihedrals of the nucleoside it can possible build with.
    # Retrieve only the keys of the angles dictionary and remove the chi dihedral as it is not of our interest here.
    backbone_list = list(json.loads(next_nucleoside.jsonObject["angles"]["dihedrals"]).keys())[:-1]

    # Get the alpha to build on, for the previous nucleoside, and then the last two to build the next nucleoside and position it correctly
    list_of_dihedrals = ["alpha" , backbone_list[-1], backbone_list[-2]]
    return list_of_dihedrals


######## FUNCTIONS USED TO GENERATE THE LIST OF THE COMPLEMENTARY NUCLEOSIDES ########
def retrieve_chemistry(chemistry : str) -> str:
    """ Parse the chemistry denominator from the nucleoside string.
        Note that the last element in the string should always be the nucleobase denominator (A, C, G, T, U).
        That is why we can safely do this. """
    ln_str = len(chemistry)
    return chemistry[:ln_str], ln_str


def retrieve_bases_list(sequence : list) -> list:
    """ retrieve the base denominator from the list of nucleic acids """
    return [x[-1] for x in sequence]


def retrieve_chemistry_list(sequence : list) -> list:
    """ retrieve the chemistry denominator from the list of nucleic acids """
    return [x[:-1] for x in sequence]


def get_complementary_bases(sequence : list, complement_dict : dict) -> list:
    """ get the complementary strand """
    return [complement_dict[x] for x in sequence]


def concatenate_chem_and_bases(chemistry : Union[str, list], bases : list) -> list:
    """ concatenate the chemistries with the bases """

    if isinstance(chemistry, list):
        nucleoside_list = [chemistry[i] + bases[i] for i in range(len(bases))]
        return nucleoside_list

    if isinstance(chemistry, str):
        nucleoside_list = [chemistry + bases[i] for i in range(len(bases))]
        return nucleoside_list


#def retrieve_base(base : str) -> str:
#    """ retrieve the base denominator for this specific nucleoside """
#    return x[-1]


def retrieve_homo_nucleosides(codex_dict_keys : list, chem_i : str, ln_str : int) -> list:
    """ Retrieves all the possible nucleobases that the chemistry can have. """

    # Outputted list of homo nucleosides
    list_of_homo_nucleosides = []

    # Slice the keys of the codex_dict_keys, meaning you only get the chemistries of all the available nucleosides
    codex_dict_keys_SLICED = [ i[:ln_str] for i in codex_dict_keys]

    # Check which index of the values that correspond to the nucleoside you are looking for
    for i, key in enumerate(codex_dict_keys_SLICED):
        if codex_dict_keys_SLICED[i] == chem_i:
            list_of_homo_nucleosides.append(codex_dict_keys[i])

    return list_of_homo_nucleosides


def assess_possible_complementary_base(chemistry : str, leadingstrand_base : str, homo_chemistry_list : list, RNA_dict : dict, DNA_dict : dict) -> str:
    """ Check which combination are valid combination with the complementary bases's dict (DNA and RNA dict) """

    dna_nucleoside = chemistry + DNA_dict[leadingstrand_base]
    rna_nucleoside = chemistry + RNA_dict[leadingstrand_base]

    # This is for A, C and G nucleosides. Since only T and U are different, we end the function here and only assess the bases with those nucleosides.
    if dna_nucleoside == rna_nucleoside:
        return dna_nucleoside

    if dna_nucleoside in homo_chemistry_list:
        return dna_nucleoside

    # If it does not exist in the list, there is only one other option, which is the rna nucleoside variant
    return rna_nucleoside


def return_chemistrycode(identifier : str) -> str:
    """ Go through the currently available files in the /Daedalus/json/ directory and read them in order to find the prompted identifier variable
        in the 'identity' list.

        Then parse the correct abbreviated chemistry identifier.
        Example ; find 'Deoxy Ribonucleic Acid' --> returns 'd'     """

    # Get the directory from the $HOME of Daedalus (where it is installed)
    from sysDaedalus import return_DAEDALUS_home
    dH = return_DAEDALUS_home()
    dirJSON = os.listdir(dH + 'json/')

    # Read all the files in the json/ directory and find the identifier that was prompted, this way we can inambiguously find the abbreviated chemical code
    for json_file in dirJSON:
        json_pathname = dH + 'json/' + json_file

        with open(json_pathname, "r") as jsonf:
            jsonContent = json.loads(json.load(jsonf)["identity"])[0]

        # if the identifier has been found, remember the name of the file we found it in
        if jsonContent == identifier:
            file_of_chemistry = json_pathname
            break


    # Iterate over the complementary codex. If you find the filename, remember the abbreviated nucleic acid code of the file 
    COMPL_CODEX = LFT3.conformations_codex
    CHECK = False
    while not CHECK:

        # Iterate over the complementary codex
        for key, value in COMPL_CODEX.items():
            if isinstance(value, list):
                for i, item in enumerate(value):
                    if value[i] == file_of_chemistry:
                        _KEY = key
                        CHECK = True

            else:
                if item == file_of_chemistry:
                    _KEY = key
                    CHECK = True

    # Return the chemistry code for the chemistry, without the nucleobase appendend, so we can add the complementary bases to it later
    return _KEY[:-1]


######## FUNCTIONS USED TO CREATE THE INPUTS FOR THE EVENTUAL DATAFRAME THAT WILL EVENTUALLY BE WRITTEN TO A PDB FORMATTED FILE ########
def return_PDB_AtomNames_or_ElementSymbol(list_of_sequence : list, identifier : str) -> list:
    """ Loads in the atom names from the json files
    The identifier is either the string "Atoms" or the string "ElementSymbol", which will parse the list of interest """

    # Initialise an empty list
    atom_list = []

    for i in range(len(list_of_sequence)):

        # If it is the first in the nucleotide sequence, meaning the last in the nucleotide sequence
        if i == 0:
            buildingblock = list_of_sequence[i]
            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            tmp_atomlist = json.loads(nucleoside["pdb_properties"][identifier])

            atom_list = atom_list + tmp_atomlist

        # If it is the last in the reversed sequence, meaning the first nucleotide
        if (i + 1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = CODEX[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            tmp_atomlist = json.loads(linker["pdb_properties"][identifier]) + json.loads(nucleoside["pdb_properties"][identifier])

            atom_list = atom_list + tmp_atomlist

            return atom_list

        if not i == 0 and not (i+1) == len(list_of_sequence):
            # since this is not the last one or the first one, just carry on as usual
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = CODEX[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            tmp_atomlist = json.loads(linker["pdb_properties"][identifier]) + json.loads(nucleoside["pdb_properties"][identifier])

            atom_list = atom_list + tmp_atomlist


def return_PDB_Sequence(list_of_sequence : list, start_of_sequence : int = 0) -> np.array :
    """ Determines the number in the sequence of the nucleotides in the strands based off on the shape of their array.
        The '+1' for the first and last nucleotide is to account for the capping of the nucleoside, here with a single hydrogen. """

    # Initialise an empty array
    sequence_array = np.array([], dtype=int)

    # Initialise a counter
    seq_count = 1 + start_of_sequence

    for i in range(len(list_of_sequence)):

        # If it is the first in the nucleotide sequence, meaning the last in the nucleotide sequence
        if i == 0:
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            nuc_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + 1
            tmp_seqarray = np.full(nuc_shape, seq_count, dtype=int)

            seq_count += 1

            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)

        # If it is the last in the reversed sequence, meaning the first nucleotide
        if (i + 1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = CODEX[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0] + 1
            tmp_seqarray = np.full(nucleotide_shape, seq_count, dtype=int)

            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)

            # Return the output of the atom_list, as this is the last nucleotide in the sequence
            return sequence_array

        # since this is not the last one or the first one, just carry on as usual
        if not i == 0 and not (i+1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = CODEX[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
            tmp_seqarray = np.full(nucleotide_shape, seq_count, dtype=int)

            seq_count += 1

            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)


def return_PDB_Residuename(list_of_sequence : list) -> list:

    # Initialise an empty array
    resname_list = []

    for i in range(len(list_of_sequence)):

        # If it is the first in the nucleotide sequence, meaning the last in the nucleotide sequence
        if i == 0:
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            nuc_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + 1
            ID = json.loads(nucleoside["identity"])[2]

            tmp_resname = [ID for i in range(nuc_shape)]
            resname_list = resname_list + tmp_resname

        # If it is the last nucleotide in the sequence
        if (i + 1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = CODEX[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0] + 1
            ID = json.loads(nucleoside["identity"])[2]

            tmp_resname = [ID for i in range(nucleotide_shape)]
            resname_list = resname_list + tmp_resname
            # Return the output of the atom_list, as this is the last nucleotide in the sequence
            return resname_list


        # since this is not the last one or the first one, just carry on as usual
        if not i == 0 and not (i+1) == len(list_of_sequence):
            buildingblock = list_of_sequence[i]

            # Make json object of nucleoside
            nucleoside = CODEX[buildingblock][0]
            with open(nucleoside, "r") as nuc:
                nucleoside = json.load(nuc)

            # Make json object of linker
            linker = CODEX[buildingblock][1]
            with open(linker, "r") as lnk:
                linker = json.load(lnk)

            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
            ID = json.loads(nucleoside["identity"])[2]

            tmp_resname = [ID for i in range(nucleotide_shape)]
            resname_list = resname_list + tmp_resname


 #### TEMPORARILY NOT IN USE
#def COMPLEMENTARY_pdb_AtomNames_or_ElementSymbol(list_of_sequence : list, identifier : str) -> list:
#    """ Loads in the atom names from the json files
#    The identifier is either the string "Atoms" or the string "ElementSymbol", which will parse the list of interest """
#
#    # Initialise an empty list
#    atom_list = []
#
#    for i in range(len(list_of_sequence)):
#
#        # If it is the first in the nucleotide sequence
#        if i == 0:
#            buildingblock = list_of_sequence[i]
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            tmp_atomlist = json.loads(nucleoside["pdb_properties"][identifier])
#
#            atom_list = atom_list + tmp_atomlist
#
#        # If it is the last in the sequence
#        if (i + 1) == len(list_of_sequence):
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            # Make json object of linker
#            linker = CODEX[buildingblock][1]
#            with open(linker, "r") as lnk:
#                linker = json.load(lnk)
#
#            tmp_atomlist = json.loads(linker["pdb_properties"][identifier]) + json.loads(nucleoside["pdb_properties"][identifier])
#
#            atom_list = atom_list + tmp_atomlist
#
#            return atom_list
#
#        if not i == 0 and not (i+1) == len(list_of_sequence):
#            # since this is not the last one or the first one, just carry on as usual
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            # Make json object of linker
#            linker = CODEX[buildingblock][1]
#            with open(linker, "r") as lnk:
#                linker = json.load(lnk)
#
#            tmp_atomlist = json.loads(linker["pdb_properties"][identifier]) + json.loads(nucleoside["pdb_properties"][identifier])
#
#            atom_list = atom_list + tmp_atomlist
#
#
#def COMPLEMENTARY_pdb_Sequence(list_of_sequence : list, start_of_sequence : int = 0) -> np.array :
#    """ Determines the number in the sequence of the nucleotides in the strands based off on the shape of their array.
#        The '+1' for the first and last nucleotide is to account for the capping of the nucleoside, here with a single hydrogen. """
#
#    # Initialise an empty array
#    sequence_array = np.array([], dtype=int)
#
#    # Initialise a counter
#    seq_count = 1 + start_of_sequence
#
#    for i in range(len(list_of_sequence)):
#
#        # If it is the first in the nucleotide sequence
#        if i == 0:
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            nuc_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + 1
#            tmp_seqarray = np.full(nuc_shape, seq_count, dtype=int)
#
#            seq_count += 1
#
#            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)
#
#        # If it is the last in the sequence
#        if (i + 1) == len(list_of_sequence):
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            # Make json object of linker
#            linker = CODEX[buildingblock][1]
#            with open(linker, "r") as lnk:
#                linker = json.load(lnk)
#
#            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0] + 1
#            tmp_seqarray = np.full(nucleotide_shape, seq_count, dtype=int)
#
#            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)
#
#            # Return the output of the atom_list, as this is the last nucleotide in the sequence
#            return sequence_array
#
#        # since this is not the last one or the first one, just carry on as usual
#        if not i == 0 and not (i+1) == len(list_of_sequence):
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            # Make json object of linker
#            linker = CODEX[buildingblock][1]
#            with open(linker, "r") as lnk:
#                linker = json.load(lnk)
#
#            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
#            tmp_seqarray = np.full(nucleotide_shape, seq_count, dtype=int)
#
#            seq_count += 1
#
#            sequence_array = np.concatenate((sequence_array, tmp_seqarray), axis=None)
#
#
#def COMPLEMENTARY_pdb_Residuename(list_of_sequence : list) -> list:
#
#    # Initialise an empty array
#    resname_list = []
#
#    for i in range(len(list_of_sequence)):
#
#        # If it is the first in the nucleotide sequence, meaning the last in the nucleotide sequence
#        if i == 0:
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            nuc_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + 1
#            ID = json.loads(nucleoside["identity"])[2]
#
#            tmp_resname = [ID for i in range(nuc_shape)]
#            resname_list = resname_list + tmp_resname
#
#        # If it is the last nucleotide in the sequence
#        if (i + 1) == len(list_of_sequence):
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            # Make json object of linker
#            linker = CODEX[buildingblock][1]
#            with open(linker, "r") as lnk:
#                linker = json.load(lnk)
#
#            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0] + 1
#            ID = json.loads(nucleoside["identity"])[2]
#
#            tmp_resname = [ID for i in range(nucleotide_shape)]
#            resname_list = resname_list + tmp_resname
#            # Return the output of the atom_list, as this is the last nucleotide in the sequence
#            return resname_list
#
#
#        # since this is not the last one or the first one, just carry on as usual
#        if not i == 0 and not (i+1) == len(list_of_sequence):
#            buildingblock = list_of_sequence[i]
#
#            # Make json object of nucleoside
#            nucleoside = CODEX[buildingblock][0]
#            with open(nucleoside, "r") as nuc:
#                nucleoside = json.load(nuc)
#
#            # Make json object of linker
#            linker = CODEX[buildingblock][1]
#            with open(linker, "r") as lnk:
#                linker = json.load(lnk)
#
#            nucleotide_shape = json.loads(nucleoside["pdb_properties"]["Shape"])[0] + json.loads(linker["pdb_properties"]["Shape"])[0]
#            ID = json.loads(nucleoside["identity"])[2]
#
#            tmp_resname = [ID for i in range(nucleotide_shape)]
#            resname_list = resname_list + tmp_resname
