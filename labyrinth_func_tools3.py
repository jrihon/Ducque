import numpy as np
import os, sys
from typing import Tuple
""" This script will function as a repository for all the atoms, angles and dihedrals that need to be parsed to build up the nucleic acid duplex"""


### WHILE THE DIHEDRAL SHOULD ALWAYS BE AROUND 180, THE ANGLES ARE OFTEN NOT A CONSTANT VALUE AND SHOULD BE INCLUDED FROM BASE-PAIR TO BASE-PAIR IN COMPLEMENT1 AND COMPLEMENT2


def retrieve_atoms_for_plane_rotation_of_complement(base1 : str, base2 : str) -> Tuple[list, list]:
    """ base1 belongs to the leading strand, base2 to the complementary base"""

    # BASE1
    if base1 == "A": atoms1 = ["N9", "C4", "C8"]
    if base1 == "C": atoms1 = ["N1", "C2", "C6"]
    if base1 == "G": atoms1 = ["N9", "C4", "C8"]
    if base1 == "T": atoms1 = ["N1", "C2", "C6"]
    if base1 == "U": atoms1 = ["N1", "C2", "C6"]

    # BASE2
    if base2 == "A": atoms2 = ["N9", "C4", "C8"]
    if base2 == "C": atoms2 = ["N1", "C2", "C6"]
    if base2 == "G": atoms2 = ["N9", "C4", "C8"]
    if base2 == "T": atoms2 = ["N1", "C2", "C6"]
    if base2 == "U": atoms2 = ["N1", "C2", "C6"]

    return atoms1, atoms2


def retrieve_atoms_for_positioning_of_complement1(base1 : str, base2 :str) -> Tuple[list, str]:
    """ base1 belongs to the leading strand, base2 to the complementary base"""

    # BASE1
    if base1 == "A": atoms1 = ["N3", "C2", "N1"]
    if base1 == "C": atoms1 = ["N1", "C2", "N3"]
    if base1 == "G": atoms1 = ["N3", "C2", "N1"]
    if base1 == "T": atoms1 = ["C2", "N3", "H3"]
    if base1 == "U": atoms1 = ["C2", "N3", "H3"]

    # BASE2
    if base2 == "A": atom2 = "N1"
    if base2 == "C": atom2 = "N3"
    if base2 == "G": atom2 = "H1"
    if base2 == "T": atom2 = "H3"
    if base2 == "U": atom2 = "H3"

    return atoms1, atom2


def retrieve_atoms_for_position_of_complement2(base1 : str, base2 :str) -> Tuple[list, str]:
    """ base1 belongs to the leading strand, base2 to the complementary base"""

    # BASE1
    if base1 == "A": atoms1 = ["C6", "N6", "H61"]
    if base1 == "C": atoms1 = ["C4", "N4", "H41"]
    if base1 == "G": atoms1 = ["N1", "C6", "O6"]
    if base1 == "T": atoms1 = ["N3", "C4", "O4"]
    if base1 == "U": atoms1 = ["N3", "C4", "O4"]

    # BASE2
    if base2 == "A": atom2 = "H61"
    if base2 == "C": atom2 = "H41"
    if base2 == "G": atom2 = "O6"
    if base2 == "T": atom2 = "O4"
    if base2 == "U": atom2 = "O4"

    return atoms1, atom2



