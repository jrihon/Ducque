import numpy as np
import os, sys, json
from typing import Tuple
""" This script will function as a repository for all the atoms, angles and dihedrals that need to be parsed to build up the nucleic acid duplex"""


codex_acidum_nucleicum = {
"dA": ["json/dna_adenosine.json", "json/dna_phosphate.json"],
"dG": ["json/dna_guanosine.json", "json/dna_phosphate.json"],
"dC": ["json/dna_cytidine.json", "json/dna_phosphate.json"],
"dT": ["json/dna_thymidine.json", "json/dna_phosphate.json"],
"rA": ["json/rna_adenosine.json", "json/rna_phosphate.json"],
"rC": ["json/rna_cytidine.json", "json/rna_phosphate.json"],
"rG": ["json/rna_guanosine.json", "json/rna_phosphate.json"],
"rU": ["json/rna_uracil.json", "json/rna_phosphate.json"],
}



complementary_codex = {
"dA": ["json/dna_adenosine_2endo.json", "json/dna_adenosine_3endo.json"],
"dC": ["json/dna_cytidine_2endo.json", "json/dna_cytidine_3endo.json"],
"dG": ["json/dna_guanosine_2endo.json", "json/dna_guanosine_3endo.json"],
"dT": ["json/dna_thymidine_2endo.json", "json/dna_thymidine_3endo.json"],
"rA": ["json/rna_adenosine_3endo.json"],
"rC": ["json/rna_cytidine_3endo.json"],
"rG": ["json/rna_guanosine_3endo.json"],
"rU": ["json/rna_uracil_3endo.json"],
}


backbone_codex = {
"DNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"RNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"Phosphate" : ["P"],
}

linker_codex = {
"Phosphate" : ["P", "OP2", "OP1"]

}

def Atom_Parsing_List(prevnuc, link, nextnuc = None) -> list:
    """ Retrieves the atoms that correspond to the correct index of the array, with which we calculate with.
    All variables are json object

    By default, nextnuc is equal to None. If no json object has been parsed into nextnuc, that means you are positioning the following linker and not the following nucleoside. """

    if nextnuc == None:
        prevChem = backbone_codex[json.loads(prevnuc.jason["identity"])[1]]
        linkChem = linker_codex[json.loads(link.jason["identity"])[0]]

        truncPrevChem = [prevChem[-3], prevChem[-2], prevChem[-1]]

        return truncPrevChem + linkChem


    if not nextnuc == None:
        prevChem = backbone_codex[json.loads(prevnuc.jason["identity"])[1]]
        linkChem = backbone_codex[json.loads(link.jason["identity"])[0]]
        nextChem = backbone_codex[json.loads(nextnuc.jason["identity"])[1]]

        if len(linkChem) == 1:
            truncPrevChem = [prevChem[-2], prevChem[-1]]

            return truncPrevChem + linkChem + nextChem

        if len(linkChem) == 2:
            truncprevChem = [prevChem[-1]]

            return truncPrevChem + linkChem + nextChem

        if len(linkChem) == 3:

            return linkChem + nextChem


def Dihedral_and_Angle_Parsing_List():
    """ Retrieves the atoms that correspond to the correct index of the array, with which we calculate with.
    All variables are json object """

    pass

def retrieve_atoms_for_plane_rotation_of_complement(base1 : str, base2 : str) -> Tuple[list, list]:
    """ base1 belongs to the leading strand, base2 to the complementary base """

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
    """ base1 belongs to the leading strand, base2 to the complementary base """

    # BASE1
    if base1 == "A": atoms1 = ["N3", "C2", "N1"]
    if base1 == "C": atoms1 = ["N1", "C2", "N3"]
    if base1 == "G": atoms1 = ["C2", "N1", "H1"]
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
    """ base1 belongs to the leading strand, base2 to the complementary base """

    # BASE1
    if base1 == "A": atoms1 = ["C6", "N6", "H61"]
    if base1 == "C": atoms1 = ["C4", "N4", "H41"]
    if base1 == "G": atoms1 = ["C5", "C6", "O6"]
    if base1 == "T": atoms1 = ["C5", "C4", "O4"]
    if base1 == "U": atoms1 = ["C5", "C4", "O4"]

    # BASE2
    if base2 == "A": atom2 = "H61"
    if base2 == "C": atom2 = "H41"
    if base2 == "G": atom2 = "O6"
    if base2 == "T": atom2 = "O4"
    if base2 == "U": atom2 = "O4"

    return atoms1, atom2


def retrieve_angles_and_dihedrals_for_initial_base_positioning(base1 : str) -> Tuple[float, float, float, float]:
    """ Only the leading strand's base is required. Normally we could hardcore in pairs, but this would disallow mismatching.

        Return in the following order :
        Q angle, Q dihedral, R angle, R dihedral            """

    Q_dihedral = 180
    R_dihedral = 180

    to_rad = (np.pi / 180)

    if base1 == "A": Q_angle, R_angle = 121.822 * to_rad, 178.802 * to_rad
    if base1 == "C": Q_angle, R_angle = 117.050 * to_rad, 176.177 * to_rad
    if base1 == "G": Q_angle, R_angle = 177.195 * to_rad, 123.466 * to_rad
    if base1 == "T": Q_angle, R_angle = 178.359 * to_rad, 122.583 * to_rad
    if base1 == "U": Q_angle, R_angle = 178.359 * to_rad, 122.583 * to_rad

    return Q_angle, Q_dihedral, R_angle, R_dihedral

