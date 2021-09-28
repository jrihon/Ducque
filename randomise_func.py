import json
import random
import sys, os
from typing import Union

import numpy as np


chemistry_dict = {
        "DNA" : ["dA", "dC", "dG", "dT"],
        "RNA" : ["rA", "rC", "rG", "rU"],
        }


def randomise_sequence(chemistry : str, length_seq : int):
    """ depends on the prompted chemistry and the length of the sequence """

    list_of_nucleotides = chemistry_dict[chemistry]

    return ", ".join(random.choice(list_of_nucleotides) for i in range(length_seq))


def randomise_chemistry(chemistry : list, sequence : list):
    """ randomise the chemistry for a given sequence """
    # Get the given sequence and chemistry into a list without any comma's
    sequence = list(map(lambda x: x.strip(","), sequence))
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    # Include the different chemistries by parsing from the dictionary, taking the first option and then cutting out the base part.. leaving only the chemistry
    list_of_chemistries = []
    for chem in chemistry:
        list_of_chemistries.append(chemistry_dict[chem][0])

    list_of_chemistries = list(map(lambda x : x[:-1], list_of_chemistries))

    # Concatenate strings by adding a random chemistry, from the prompted list, 
    tmpseq = []
    for nucleotide in sequence:
        nucl = random.choice(list_of_chemistries) + nucleotide
        tmpseq.append(nucl)

    return ", ".join(tmpseq)


def randomise_sequence_and_chemistry(chemistry : list, length_seq : int):
    """ randomise both the sequence and for a set of given chemistry """
    # Get the given chemistry without the comma's in the strings
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    list_of_nucleotides = []
    for chem in chemistry:
        list_of_nucleotides.extend(chemistry_dict[chem])

    return ", ".join(random.choice(list_of_nucleotides) for i in range(length_seq))


def write_out_complementary_sequence(compl_seq : Union[str, list]) -> str:

    if isinstance(compl_seq):
        compl_seq = list(map(lambda x: x.strip()))
        output_sequence = ", ".join(compl_seq)
        return output_sequence

    return compl_seq
