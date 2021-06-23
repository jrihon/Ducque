import random
import sys, os
import json
import numpy as np
from typing import Union

""" This script is used to randomise an outputted sequence based on a given input of parameters. """

chemistry_dict = {
        "DNA" : ["dA", "dC", "dG", "dT"],
        "RNA" : ["rA", "rC", "rG", "rU"],
        }


def randomise_sequence(chemistry : str, length_seq : int):
    """ depends on the prompted chemistry and the length of the sequence """

    list_of_nucleotides = chemistry_dict[chemistry]

    output_sequence = ", ".join(random.choice(list_of_nucleotides) for i in range(length_seq))

    with open("./sequence_testing.in", "w") as seq:
        seq.write(output_sequence + "\n")


def randomise_chemistry(chemistry : list, sequence : list):
    """ randomise the chemistry for a given sequence """
    # Get the given sequence and chemistry into a list without any comma's
    sequence = list(map(lambda x: x.strip(","), sequence))
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    # Include the different chemistries by parsing from the dictionary and then cutting out the base part
    list_of_chemistries = []
    for chem in chemistry:
        list_of_chemistries.append(chemistry_dict[chem][0])

    list_of_chemistries = list(map(lambda x : x[:-1], list_of_chemistries))

    # Concatenate strings by adding a random chemistry, from the prompted list, 
    tmpseq = []
    for nucleotide in sequence:
        nucl = random.choice(list_of_chemistries) + nucleotide
        tmpseq.append(nucl)

    output_sequence = ", ".join(tmpseq)

    with open("./sequence_testing.in", "w") as seq:
        seq.write(output_sequence + "\n")



def randomise_sequence_and_chemistry(chemistry : list, length_seq : int):
    """ randomise both the sequence and for a set of given chemistry """
    # Get the given chemistry without the comma's in the strings
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    list_of_nucleotides = []
    for chem in chemistry:
        list_of_nucleotides.extend(chemistry_dict[chem])

    output_sequence = ", ".join(random.choice(list_of_nucleotides) for i in range(length_seq))

    with open("./sequence_testing.in", "w") as seq:
        seq.write(output_sequence + "\n")


def randomiser(chemistry : Union[str, list], length_seq : int, sequence : list):
    """ the 'MAIN' function in the Daedalus.randomise script """

    # Could do with a try: these instances -> Except raise error. mutuallyexclusive.
    if length_seq != 0 and sequence is not None:
        print("Both a sequence and its length are specified, but these are mutually exclusive. Please check the given input file.\n")
        sys.exit(0)

    # if randomise sequence, so if there is only one chemistry given
    if isinstance(chemistry, str) and length_seq:
        return randomise_sequence(chemistry, length_seq)

    # if randomise chemistry for a given sequence, so a sequence must be adhered to, but not the chemistry.
    # can choose to input a set of chemistry or all chemistry
    if chemistry is not None and sequence is not None:
        return randomise_chemistry(chemistry, sequence)

    # if randomise both the sequence and the chemistry of the sequence
    # justfuckmeupfam
    if chemistry is not None and length_seq > 0:
        return randomise_sequence_and_chemistry(chemistry, length_seq)

