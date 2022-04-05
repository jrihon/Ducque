import json
import random
import sys, os
from typing import Union

import numpy as np
import fundaments
from labyrinth_func_tools2 import return_chemistrycode
from labyrinth_repository import codex_acidum_nucleicum

CODEX = codex_acidum_nucleicum

chemistry_dict = {
        "DNA" : ["dA", "dC", "dG", "dT"],
        "RNA" : ["rA", "rC", "rG", "rU"],
        }

def join_chemistry_with_sequence(chemistry : str, sequence) -> str:
    """ For the lazy people, concatenate the chemistry with the sequence they want. Only works for single prompted chemistries"""

    sequence = list(map(lambda x: x.strip(","), sequence))

    NUC_ID = fundaments.check_if_chemistry_is_valid(chemistry)

    # Get the abbreviated chemistry code
    abbrCode = return_chemistrycode(NUC_ID)

    return ", ".join(abbrCode + Base for Base in sequence)


def randomise_sequence(chemistry : str, length_seq : int) -> str:
    """ depends on the prompted chemistry and the length of the sequence """

    NUC_ID = fundaments.check_if_chemistry_is_valid(chemistry)
    abbrCode = return_chemistrycode(NUC_ID)

    all_abbrCodes = list(CODEX.keys())

    list_of_possible_nucleotides = []
    for i, code in enumerate(all_abbrCodes):
        if code[:-1] == abbrCode:
            list_of_possible_nucleotides.append(all_abbrCodes[i])

    return ", ".join(random.choice(list_of_possible_nucleotides) for i in range(length_seq))


def randomise_chemistry(chemistry : list, sequence : list) -> str:
    """ randomise the chemistry for a given sequence """

    # Get the given sequence and chemistry into a list without any comma's
    sequence = list(map(lambda x: x.strip(","), sequence))
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    # Include the different chemistries by parsing from the dictionary, taking the first option and then cutting out the base part. leaving only the chemistry
    list_of_chemistries = []
    for chem in chemistry:
        NUC_ID = fundaments.check_if_chemistry_is_valid(chem)
        abbrCode = return_chemistrycode(NUC_ID)
        list_of_chemistries.append(abbrCode)


    # Concatenate strings by adding a random chemistry, from the prompted list, 
    tmpseq = []
    for nucleotide in sequence:
        nucl = random.choice(list_of_chemistries) + nucleotide
        tmpseq.append(nucl)

    return ", ".join(tmpseq)


def randomise_sequence_and_chemistry(chemistry : list, length_seq : int) -> str:
    """ randomise both the sequence and for a set of given chemistry """
    # Get the given chemistry without the comma's in the strings
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    # Get all the abbreviated chemistry codes from the prompted chemistries you want to randomise
    list_of_chemistries = []
    for chem in chemistry:
        NUC_ID = fundaments.check_if_chemistry_is_valid(chem)
        abbrCode = return_chemistrycode(NUC_ID)
        list_of_chemistries.append(abbrCode)

    all_abbrCodes = list(CODEX.keys())

    list_of_possible_nucleotides = []
    for i, code in enumerate(all_abbrCodes):
        slicedCode = code[:-1]

        if slicedCode in list_of_chemistries:
            list_of_possible_nucleotides.append(all_abbrCodes[i])


    return ", ".join(random.choice(list_of_possible_nucleotides) for i in range(length_seq))


def write_out_complementary_sequence(compl_seq : Union[str, list]) -> str:
    """ If it is a list, parse the list correctly and output it correctly.
        If it is just a string, output it as as it is now. """

    if isinstance(compl_seq, list):
        compl_seq = list(map(lambda x: x.strip(","), compl_seq))
        output_sequence = ", ".join(compl_seq)
        return output_sequence

    return compl_seq
