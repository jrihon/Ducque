import sys, os
import json
import numpy as np
from typing import Union

import randomise_func as RF


""" This script is used to randomise an outputted sequence based on a given input of parameters. """

def write_sequence_to_file(outFile : str, outputSequence : str, complement : str):

    fname = outFile + ".random_in"
    with open("./" + fname, "w") as seq:
        seq.write(f"--sequence {outputSequence}\n--complement {complement}\n")

    print(f"Writing to {fname} . \n\n")


def randomiser(chemistry : Union[str, list], lengthSequence : int, sequence : list, complementSequence : list, outFile : str):
    """ the 'MAIN' function in the Daedalus.randomise script """

    # RANDOMISE A SEQUENCE FOR A SINGLE CHEMISTRY
    if isinstance(chemistry, str) and lengthSequence:
        outputSequence = RF.randomise_sequence(chemistry, lengthSequence)
        complementSequence = RF.write_out_complementary_sequence(complementSequence)

        write_sequence_to_file(outFile, outputSequence, complementSequence)
        return

    # IF GIVEN A SINGLE CHEMISTRY AND A SEQUENCE, CONCATENATE THE GIVEN CHEMISTRY AND THE GIVEN SEQUENCE.
    # NB : this is not randomisation, just a bit of laziness
    if isinstance(chemistry, str) and sequence:
        outputSequence = RF.join_chemistry_with_sequence(chemistry, sequence)
        complementSequence = RF.write_out_complementary_sequence(complementSequence)

        write_sequence_to_file(outFile, outputSequence, complementSequence)
        return

    # RANDOMISE A SEQUENCE'S CHEMISTRIES FROM A GIVEN SEQUENCE AND A GIVEN LIST OF CHEMISTRIES
    if chemistry is not None and sequence is not None:
        outputSequence = RF.randomise_chemistry(chemistry, sequence)
        complementSequence = RF.write_out_complementary_sequence(complementSequence)

        write_sequence_to_file(outFile, outputSequence, complementSequence)
        return

    # RANDOMISE BOTH THE CHEMISTRIES AND THE SEQUENCE FOR A GIVEN LENGTH AND A GIVEN LIST OF CHEMISTRIES
    if chemistry is not None and lengthSequence > 0:
        outputSequence = RF.randomise_sequence_and_chemistry(chemistry, lengthSequence)
        complementSequence = RF.write_out_complementary_sequence(complementSequence)

        write_sequence_to_file(outFile, outputSequence, complementSequence)
        return

