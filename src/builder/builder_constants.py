from typing import Tuple, List
from numpy import pi



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


def retrieve_atoms_for_positioning_of_complement1(base1 : str, base2 :str) -> Tuple[List[str], str, float]:
    """ base1 belongs to the leading strand, base2 to the complementary base """

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


def retrieve_atoms_for_position_of_complement2(base1 : str, base2 :str) -> Tuple[List[str], str, float]:
    """ base1 belongs to the leading strand, base2 to the complementary base """

    if base1 == "A": atoms1, distance = ["C5", "C6", "N1"], 3.749
    if base1 == "C": atoms1, distance = ["C5", "C4", "N3"], 3.698
    if base1 == "G": atoms1, distance = ["C5", "C6", "N1"], 3.745
    if base1 == "T": atoms1, distance = ["C5", "C4", "N3"], 3.745
    if base1 == "U": atoms1, distance = ["C5", "C4", "N3"], 3.756

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

    to_rad = (pi / 180)

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
