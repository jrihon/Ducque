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

    # if the chemistry is a string, there is only one chemistry given.
    if isinstance(chemistry, str) and length_seq:
        return randomise_sequence(chemistry, length_seq)

    # if given a single chemistry and a sequence, do nothing. If you already write out the specific sequence, might as well add the chemistry
    if isinstance(chemistry, str) and sequence:
        print("\nThat's pretty lazy bro. You've given a sequence, but can't add the single chemistry denotator in front of it? bruh.\n")
        return

    # if randomise chemistry for a given sequence, we adhere to a given sequence and can randomise with the given chemistries
    if chemistry is not None and sequence is not None:
        return randomise_chemistry(chemistry, sequence)

    # if you've reach this part, you can randomise the sequence and the chemistry (with the given inputs)
    if chemistry is not None and length_seq > 0:
        return randomise_sequence_and_chemistry(chemistry, length_seq)

