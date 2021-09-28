import sys, os
import json
import numpy as np
from typing import Union

import randomise_func as RF
""" This script is used to randomise an outputted sequence based on a given input of parameters. """


def write_sequence_to_file(output_sequence : str, complement : str):

    fname = "sequence_testing.in"
    with open("./" + fname, "w") as seq:
        seq.write("--sequence " + output_sequence + "\n--complement " + complement + "\n")

    print("Writing to " + fname + ". \n\n")


def randomiser(chemistry : Union[str, list], length_seq : int, sequence : list, compl_seq : list):
    """ the 'MAIN' function in the Daedalus.randomise script """

    # if the chemistry is a string, there is only one chemistry given.
    if isinstance(chemistry, str) and length_seq:
        output_sequence = RF.randomise_sequence(chemistry, length_seq)
        compl_seq = RF.write_out_complementary_sequence(compl_seq)

        write_sequence_to_file(output_sequence, compl_seq)
        return

    # if given a single chemistry and a sequence, do nothing. If you already write out the specific sequence, might as well add the chemistry
    if isinstance(chemistry, str) and sequence:
        print("\nThat's pretty lazy bro. You've given a sequence, but can't add the single chemistry denotator in front of it? bruh.\n")
        return

    # if randomise chemistry for a given sequence, we adhere to a given sequence and can randomise with the given chemistries
    if chemistry is not None and sequence is not None:
        output_sequence = RF.randomise_chemistry(chemistry, sequence)
        compl_seq = RF.write_out_complementary_sequence(compl_seq)

        write_sequence_to_file(output_sequence, compl_seq)
        return

    # if you've reach this part, you can randomise the sequence and the chemistry (with the given inputs)
    if chemistry is not None and length_seq > 0:
        output_sequence = RF.randomise_sequence_and_chemistry(chemistry, length_seq)
        compl_seq = RF.write_out_complementary_sequence(compl_seq)

        write_sequence_to_file(output_sequence, compl_seq)
        return

