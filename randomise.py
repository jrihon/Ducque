import sys, os
import json
import numpy as np
from typing import Union

import randomise_func as RF


""" This script is used to randomise an outputted sequence based on a given input of parameters. """

def write_sequence_to_file(output_sequence : str, complement : str):

    fname = "sequence_testing.in"
    with open("./" + fname, "w") as seq:
        seq.write(f"--sequence {output_sequence}\n--complement {complement}\n")

    print(f"Writing to {fname}. \n\n")


def randomiser(chemistry : Union[str, list], length_seq : int, sequence : list, compl_seq : list):
    """ the 'MAIN' function in the Daedalus.randomise script """

    # RANDOMISE A SEQUENCE FOR A SINGLE CHEMISTRY
    if isinstance(chemistry, str) and length_seq:
        output_sequence = RF.randomise_sequence(chemistry, length_seq)
        compl_seq = RF.write_out_complementary_sequence(compl_seq)

        write_sequence_to_file(output_sequence, compl_seq)
        return

    # IF GIVEN A SINGLE CHEMISTRY AND A SEQUENCE, CONCATENATE THE GIVEN CHEMISTRY AND THE GIVEN SEQUENCE.
    # NB : this is not randomisation, just a bit of laziness
    if isinstance(chemistry, str) and sequence:
        output_sequence = RF.join_chemistry_with_sequence(chemistry, sequence)
        compl_seq = RF.write_out_complementary_sequence(compl_seq)

        write_sequence_to_file(output_sequence, compl_seq)
        return

    # RANDOMISE A SEQUENCE'S CHEMISTRIES FROM A GIVEN SEQUENCE AND A GIVEN LIST OF CHEMISTRIES
    if chemistry is not None and sequence is not None:
        output_sequence = RF.randomise_chemistry(chemistry, sequence)
        compl_seq = RF.write_out_complementary_sequence(compl_seq)

        write_sequence_to_file(output_sequence, compl_seq)
        return

    # RANDOMISE BOTH THE CHEMISTRIES AND THE SEQUENCE FOR A GIVEN LENGTH AND A GIVEN LIST OF CHEMISTRIES
    if chemistry is not None and length_seq > 0:
        output_sequence = RF.randomise_sequence_and_chemistry(chemistry, length_seq)
        compl_seq = RF.write_out_complementary_sequence(compl_seq)

        write_sequence_to_file(output_sequence, compl_seq)
        return

