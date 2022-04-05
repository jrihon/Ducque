import numpy as np
import os, sys, json
from typing import Tuple

import sysDaedalus
""" This script will function as a repository for all the atoms, angles and dihedrals that need to be parsed to build up the nucleic acid duplex"""

dh = sysDaedalus.return_DAEDALUS_home()

codex_acidum_nucleicum = {
"dA" : [ dh + "json/dna_adenosine_2endo.json", dh + "json/dna_phosphate.json"],
"dG" : [ dh + "json/dna_guanosine_2endo.json", dh + "json/dna_phosphate.json"],
"dC" : [ dh + "json/dna_cytidine_2endo.json", dh + "json/dna_phosphate.json"],
"dT" : [ dh + "json/dna_thymidine_2endo.json", dh + "json/dna_phosphate.json"],
"rA" : [ dh + "json/rna_adenosine_3endo.json", dh + "json/rna_phosphate.json"],
"rC" : [ dh + "json/rna_cytidine_3endo.json", dh + "json/rna_phosphate.json"],
"rG" : [ dh + "json/rna_guanosine_3endo.json", dh + "json/rna_phosphate.json"],
"rU" : [ dh + "json/rna_uracil_3endo.json", dh + "json/rna_phosphate.json"],
"ddA" : [ dh + "json/b-homodna_adenosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"ddG" : [ dh + "json/b-homodna_guanosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"ddC" : [ dh + "json/b-homodna_cytidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"ddT" : [ dh + "json/b-homodna_thymidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hA" : [ dh + "json/hna_adenosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hG" : [ dh + "json/hna_guanosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hC" : [ dh + "json/hna_cytidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hT" : [ dh + "json/hna_thymidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"xA" : [ dh + "json/xylo_adenosine_3endo.json", dh + "json/dna_phosphate.json"],
"xG" : [ dh + "json/xylo_guanosine_3endo.json", dh + "json/dna_phosphate.json"],
"xC" : [ dh + "json/xylo_cytidine_3endo.json", dh + "json/dna_phosphate.json"],
"xU" : [ dh + "json/xylo_uracil_3endo.json", dh + "json/dna_phosphate.json"],
"2MA" : [ dh + "json/2-ome-rna_adenosine_3endo.json", dh + "json/dna_phosphate.json"],
"2MG" : [ dh + "json/2-ome-rna_guanosine_3endo.json", dh + "json/dna_phosphate.json"],
"2MC" : [ dh + "json/2-ome-rna_cytidine_3endo.json", dh + "json/dna_phosphate.json"],
"2MU" : [ dh + "json/2-ome-rna_uracil_3endo.json", dh + "json/dna_phosphate.json"],
"cA" : [ dh + "json/cena_adenosine_3endo.json", dh + "json/rna_phosphate.json"],
"cG" : [ dh + "json/cena_guanosine_3endo.json", dh + "json/rna_phosphate.json"],
"cC" : [ dh + "json/cena_cytidine_3endo.json", dh + "json/rna_phosphate.json"],
"cT" : [ dh + "json/cena_thymidine_3endo.json", dh + "json/rna_phosphate.json"],
"dxA" : [dh + "json/dxylo_adenosine_3endo.json", dh + "json/dna_phosphate.json"],
"dxG" : [dh + "json/dxylo_guanosine_3endo.json", dh + "json/dna_phosphate.json"],
"dxC" : [dh + "json/dxylo_cytidine_3endo.json", dh + "json/dna_phosphate.json"],
"dxT" : [dh + "json/dxylo_thymidine_3endo.json", dh + "json/dna_phosphate.json"],
"mA" : [dh + "json/mna_adenosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"mG" : [dh + "json/mna_guanosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"mC" : [dh + "json/mna_cytidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"mT" : [dh + "json/mna_thymidine_1-4chair.json", dh + "json/dna_phosphate.json"],
}



conformations_codex = {
"dA": [ dh + "json/dna_adenosine_2endo.json", dh + "json/dna_adenosine_3endo.json"],
"dC": [ dh + "json/dna_cytidine_2endo.json", dh + "json/dna_cytidine_3endo.json"],
"dG": [ dh + "json/dna_guanosine_2endo.json", dh + "json/dna_guanosine_3endo.json"],
"dT": [ dh + "json/dna_thymidine_2endo.json", dh + "json/dna_thymidine_3endo.json"],
"rA": [ dh + "json/rna_adenosine_3endo.json"],
"rC": [ dh + "json/rna_cytidine_3endo.json"],
"rG": [ dh + "json/rna_guanosine_3endo.json"],
"rU": [ dh + "json/rna_uracil_3endo.json"],
"ddA" : [ dh + "json/b-homodna_adenosine_1-4chair.json"], # "json/b-homodna_adenosine_4-1chair.json"],
"ddG" : [ dh + "json/b-homodna_guanosine_1-4chair.json"], #"json/b-homodna_guanosine_4-1chair.json"],
"ddC" : [ dh + "json/b-homodna_cytidine_1-4chair.json"], #"json/b-homodna_cytosine_4-1chair.json"],
"ddT" : [ dh + "json/b-homodna_thymidine_1-4chair.json"], # "json/b-homodna_thymidine_4-1chair.json",],
"hA" : [ dh + "json/hna_adenosine_1-4chair.json"], # "json/hna_adenosine_4-1chair.json"],
"hG" : [ dh + "json/hna_guanosine_1-4chair.json"], # "json/hna_guanosine_4-1chair.json"],
"hC" : [ dh + "json/hna_cytidine_1-4chair.json"], # "json/hna_cytidine_4-1chair.json"],
"hT" : [ dh + "json/hna_thymidine_1-4chair.json"], # "json/hna_thymidine_4-1chair.json"],
"xA" : [ dh + "json/xylo_adenosine_3endo.json"],
"xG" : [ dh + "json/xylo_guanosine_3endo.json"],
"xC" : [ dh + "json/xylo_cytidine_3endo.json"],
"xU" : [ dh + "json/xylo_uracil_3endo.json"],
"2MA" : [ dh + "json/2-ome-rna_adenosine_3endo.json"],
"2MG" : [ dh + "json/2-ome-rna_guanosine_3endo.json"],
"2MC" : [ dh + "json/2-ome-rna_cytidine_3endo.json"],
"2MU" : [ dh + "json/2-ome-rna_uracil_3endo.json"],
"cA": [ dh + "json/cena_adenosine_3endo.json"],
"cG": [ dh + "json/cena_guanosine_3endo.json"],
"cC": [ dh + "json/cena_cytidine_3endo.json"],
"cT": [ dh + "json/cena_thymidine_3endo.json"],
"dxA" : [dh + "json/dxylo_adenosine_3endo.json"],
"dxG" : [dh + "json/dxylo_guanosine_3endo.json"],
"dxC" : [dh + "json/dxylo_cytidine_3endo.json"],
"dxT" : [dh + "json/dxylo_thymidine_3endo.json"],
"mA" : [dh + "json/mna_adenosine_1-4chair.json"],
"mG" : [dh + "json/mna_guanosine_1-4chair.json"],
"mC" : [dh + "json/mna_cytidine_1-4chair.json"],
"mT" : [dh + "json/mna_thymidine_1-4chair.json"],
}


backbone_codex = {
"Phosphate" : ["P"],
"DNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"RNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"b-homoDNA" : ["O4'", "C4'", "C5'", "C6'", "O6'"],
"HNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"Xylo" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"2-OMe-RNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"CeNA" : ["O3'", "C3'", "C4'", "C7'", "O7'"],
"dXylo" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"MNA" : ["N3'", "C4'", "C5'", "C6'", "O6'"],

                }

linker_codex = {
"Phosphate" : ["P", "OP2", "OP1"],

                }





def Atom_Parsing_List(prevnuc, link, nextnuc = None) -> list:
    """ Retrieves the atoms that correspond to the correct index of the array, with which we calculate with.
    All variables are json object

    By default, nextnuc is equal to None. If no json object has been parsed into nextnuc, that means you are positioning the following linker and not the following nucleoside. """

    # For when adding the linker moiety to the nucleoside
    if nextnuc == None:
        prevChem = backbone_codex[json.loads(prevnuc.jsonObject["identity"])[1]]
        linkChem = linker_codex[json.loads(link.jsonObject["identity"])[0]]

        truncPrevChem = [prevChem[-3], prevChem[-2], prevChem[-1]]

        return truncPrevChem + linkChem


    # For when building a subsequent nucleoside to the current nucleotide
    if not nextnuc == None:
        prevChem = backbone_codex[json.loads(prevnuc.jsonObject["identity"])[1]]
        linkChem = backbone_codex[json.loads(link.jsonObject["identity"])[0]]
        nextChem = backbone_codex[json.loads(nextnuc.jsonObject["identity"])[1]]

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

    # Set distance between the two bases (1)
#    First try
#    distance = 1.81
#
#    # BASE1
#    if base1 == "A": atoms1 = ["N3", "C2", "N1"]
#    if base1 == "C": atoms1 = ["N1", "C2", "N3"]
#    if base1 == "G": atoms1 = ["C2", "N1", "H1"]
#    if base1 == "T": atoms1 = ["C2", "N3", "H3"]
#    if base1 == "U": atoms1 = ["C2", "N3", "H3"]
#
#    # BASE2
#    if base2 == "A": atom2 = "N1"
#    if base2 == "C": atom2 = "N3"
#    if base2 == "G": atom2 = "H1"
#    if base2 == "T": atom2 = "H3"
#    if base2 == "U": atom2 = "H3"
#
#    return atoms1, atom2, distance
#
#

#    Second and third try
    # BASE1
    if base1 == "A": atoms1, distance = ["N3", "C2", "N1"], 2.90
    if base1 == "C": atoms1, distance = ["N1", "C2", "N3"], 2.87
    if base1 == "G": atoms1, distance = ["N3", "C2", "N1"], 2.87
    if base1 == "T": atoms1, distance = ["N1", "C2", "N3"], 2.90
    if base1 == "U": atoms1, distance = ["N1", "C2", "N3"], 2.90

    # BASE2
    if base2 == "A": atom2 = "N1"
    if base2 == "C": atom2 = "N3"
    if base2 == "G": atom2 = "N1"
    if base2 == "T": atom2 = "N3"
    if base2 == "U": atom2 = "N3"

    return atoms1, atom2, distance


def retrieve_atoms_for_position_of_complement2(base1 : str, base2 :str) -> Tuple[list, str]:
    """ base1 belongs to the leading strand, base2 to the complementary base """
    # Set distance between the two bases (2)
#    distance = 1.87
#    distance = 2.82

    # Do not target hydrogens for base pairing, as with optimized orbitals, the amine groups on the nucleobase are not planar any more and that makes it more difficult.

    # BASE1
#    First try
#    if base1 == "A": atoms1 = ["C6", "N6", "H61"]
#    if base1 == "C": atoms1 = ["C4", "N4", "H41"]
#    Second try
#    if base1 == "A": atoms1 = ["C5", "C6", "N6"]
#    if base1 == "C": atoms1 = ["C5", "C4", "N4"]
#    if base1 == "G": atoms1 = ["C5", "C6", "O6"]
#    if base1 == "T": atoms1 = ["C5", "C4", "O4"]
#    if base1 == "U": atoms1 = ["C5", "C4", "O4"]
#    Third Try
    if base1 == "A": atoms1, distance = ["C5", "C6", "N1"], 3.749
    if base1 == "C": atoms1, distance = ["C5", "C4", "N3"], 3.698
    if base1 == "G": atoms1, distance = ["C5", "C6", "N1"], 3.745
    if base1 == "T": atoms1, distance = ["C5", "C4", "N3"], 3.745
    if base1 == "U": atoms1, distance = ["C5", "C4", "N3"], 3.756

    # BASE2
#    First try
#    if base2 == "A": atom2 = "H61"
#    if base2 == "C": atom2 = "H41"
#    Second try
#    if base2 == "A": atom2 = "N6"
#    if base2 == "C": atom2 = "N4"
#    if base2 == "G": atom2 = "O6"
#    if base2 == "T": atom2 = "O4"
#    if base2 == "U": atom2 = "O4"
#    Third Try
    if base2 == "A": atom2 = "C6"
    if base2 == "C": atom2 = "C4"
    if base2 == "G": atom2 = "C6"
    if base2 == "T": atom2 = "C4"
    if base2 == "U": atom2 = "C4"

    return atoms1, atom2, distance


def retrieve_angles_and_dihedrals_for_initial_base_positioning(nucleobase : str) -> Tuple[float, float, float, float]:
    """ Only the leading strand's base is required. Normally we could hardcore in pairs, but this would disallow mismatching.

        Return in the following order :
        Q angle, Q dihedral, R angle, R dihedral            """

    # Regular Watson Crick Franklin base pairing
    Q_dihedral = 180
    R_dihedral = 180

    to_rad = (np.pi / 180)

    # First try
#    if base1 == "A": Q_angle, R_angle = 121.822 * to_rad, 178.802 * to_rad
#    if base1 == "G": Q_angle, R_angle = 177.195 * to_rad, 123.466 * to_rad
#    if base1 == "C": Q_angle, R_angle = 117.050 * to_rad, 176.177 * to_rad
#    if base1 == "T": Q_angle, R_angle = 178.359 * to_rad, 122.583 * to_rad
#    if base1 == "U": Q_angle, R_angle = 178.359 * to_rad, 122.583 * to_rad

    # Second try
#    if base1 == "A": Q_angle, R_angle = 121.822 * to_rad, 119.999 * to_rad
#    if base1 == "G": Q_angle, R_angle = 119.185 * to_rad, 123.466 * to_rad
#    if base1 == "C": Q_angle, R_angle = 117.050 * to_rad, 120.123 * to_rad
#    if base1 == "T": Q_angle, R_angle = 115.789 * to_rad, 122.583 * to_rad
#    if base1 == "U": Q_angle, R_angle = 115.789 * to_rad, 122.583 * to_rad

    # Third try
    if nucleobase == "A": Q_angle, R_angle = 121.822 * to_rad, 100.942 * to_rad
    if nucleobase == "C": Q_angle, R_angle = 117.050 * to_rad, 101.381 * to_rad
    if nucleobase == "G": Q_angle, R_angle = 119.185 * to_rad,  97.809 * to_rad
    if nucleobase == "T": Q_angle, R_angle = 115.789 * to_rad,  99.786 * to_rad
    if nucleobase == "U": Q_angle, R_angle = 115.789 * to_rad,  99.786 * to_rad
    return Q_angle, Q_dihedral, R_angle, R_dihedral


def retrieve_atom_for_direction_axis(nucleobase : str) -> str:
    """ To create a straight axis for turning when using tilting the axis to get a better fit """
    if nucleobase == "A": return "C4"
    if nucleobase == "C": return "C6"
    if nucleobase == "G": return "C4"
    if nucleobase == "T": return "C6"
    if nucleobase == "U": return "C6"
