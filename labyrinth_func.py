import numpy as np
import pandas as pd
import json, sys

import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2


"""
Labyrinth_func.py
The script that contains the classes and all the functions that concatenate the workflow of consecutively adding the linker and nucleotides.
"""

# CLASSES
class Nucleoside:

    def __init__(self, jsonfile):

        with open(jsonfile, "r") as jsonf:
            self.jason = json.load(jsonf)

        self.splitted = jsonfile.split(".")[0]
        self.array =  np.asarray(json.loads(self.jason["pdb_properties"]["Coordinates"]), dtype=float)
        self.atom_list = json.loads(self.jason["pdb_properties"]["Atoms"])


    def get_alpha(self) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["alpha"])

    def get_beta(self) -> float:
        # because beta is still inside a the backbone dictionary, we need to load it again 
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["beta"])

    def get_zeta(self) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["zeta"])

    def get_epsilon(self) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["epsilon"])


class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_COP(self) -> float:
        # since the angles are inside the angles dictionary, we don"t need to load string again
        return float(self.jason["Angles"]["C5_O5_P"]) * (np.pi / 180)

    def get_shape(self) -> int:
        # returns the size of the linker shape, but only the first value
        return int(json.loads(self.jason["pdb_properties"]["Shape"])[0])

    def get_OPO2(self) -> float:
        return float(self.jason["Angles"]["O5_P_OP2"]) * (np.pi / 180)

    def get_OPO1(self) -> float:
        return float(self.jason["Angles"]["O5_P_OP1"]) * (np.pi / 180)

    def get_O3PO5(self) -> float:
        return float(self.jason["Angles"]["O5_P_O3"]) * (np.pi / 180)

    def get_OP2_dihedral(self) -> float:
        return float(self.jason["Dihedrals"]["dihedral_oxygen_OP2"])

    def get_OP1_dihedral(self) -> float:
        return float(self.jason["Dihedrals"]["dihedral_oxygen_OP1"])


## FUNCTIONS THAT ARE MEANT TO BYPASS THE ITERATIVE CODING AND RESULT IN ONLY THE NECESSARY RESULTS IN THE LABYRINTH.PY : ARCHITECTURE()
def generate_vector_of_interest(angle : float, dihedral : float, atom_array : np.array) -> np.array:
    """ This function generates a single vector for which there exists only one angle and dihedral.
    The atom_array contains the three first atoms in the sequence that make up the dihedral.
    Example = if the sequence of a dihedral is C4' - C5' - O5' - P (beta backbone), then the atom_array is [O5', C5', C4']
    """
    #print(1)
    # Get the vector to rotate the cone vector onto. This is done on the middle two atoms of the sequence
    # Example : [O5'] minus (-) [C5'] results in a vector that goes [C5' -> O5']
    p0 = LFT1.return_normalized(atom_array[0] - atom_array[1])

    # Generate the cone vector
    cone_vector = LFT1.generate_cone_vector(angle)

    # Get the quaternion that corresponds to the desired rotation<
    # Note, if 'vector_to_from" is empty, it defaults to the Z-axis, around which the cone is generated np.array([0,0,1])
    quaternion = LFT1.get_quaternion(p0 * -1.0)

    # Rotate the cone vector
    rotated_cone_vector = LFT1.rotate_with_quaternion(quaternion, cone_vector)

    # Calculate the possible dihedrals we have with the first three atoms and the cone vector. Remember: the cone vector represents the position of the fourth atom
    range_of_dihedrals = LFT1.dihedral_array(atom_array, rotated_cone_vector)

    # Now extrapolate in between the range_of_dihedrals to get the correct theta angle to generate the single vector of interest.
    # The theta angle corresponds to the angle (spherical coordinates r, theta, phi) at which we need to generate our vector at.
    theta_interpolate = LFT1.get_interpolated_dihedral(range_of_dihedrals, dihedral)

    # Generate the correct vector of interest
    single_vector = LFT1.generate_and_rotate_single_vector_QUAT(theta_interpolate, angle, quaternion)

    return single_vector

def position_phosphate_linker(nucleoside, linker) -> np.ndarray:
    """
    This function positions and rotates the phosphate linker onto the nucleoside of interest.
    First the location of the phosphorus atom is calculated for.
    Secondly, the first rotation is carried out based on the P -> OP2 vector.
    Lastly, the second rotation is carried out based on the P -> OP1 vector.

    nucleoside is a json object
    linker is a json object
    """
    # Retrieve the vectors of the atoms that make up the dihedral you research
    # Dihedral C4' - C5' - O5' - P
    id_O5 = LFT2.retrieve_atom_index(nucleoside, "O5'")
    v2 = nucleoside.array[id_O5]

    id_C5 = LFT2.retrieve_atom_index(nucleoside, "C5'")
    v1 = nucleoside.array[id_C5]

    id_C4 = LFT2.retrieve_atom_index(nucleoside, "C4'")
    v0 = nucleoside.array[id_C4]

    # Find the vector that corresponds to the O5' -> P vector
    single_vector1 = generate_vector_of_interest(linker.get_COP(), nucleoside.get_beta(), [v2, v1, v0])

    # Add the single_vector1 to O5' atom, to denote the location of P. Then move the linker array to that location
    link = LFT1.position_phosphate(v2, single_vector1, linker.array)   # the linker is in the correct position  

    ## We will rotate the linker twice, to get the correct orientation of the linker in 3D space
    # Dihedral C5' - O5' - P - OP2
    id_P = LFT2.retrieve_atom_index(linker, "P")
    v3 = link[id_P]

    id_OP2 = LFT2.retrieve_atom_index(linker, "OP2")
    v4 = link[id_OP2]

    single_vector2 = generate_vector_of_interest(linker.get_OPO2(), linker.get_OP2_dihedral(), [v3, v2, v1])

    # This is the distance from the positioned phosphorus atom to the origin
    p_to_origin = v3

    # move the linker to the origin, by positioning the phosphorus at [0,0,0]
    link_to_origin = LFT1.move_vector_to_origin(link, p_to_origin)

    # Rotate the linker a first time
    # Define the vector that goes from P to OP2 and normalize it
    P_OP2 = LFT1.return_normalized(v4 - v3)

    # Get quaternion to rotate the linker a first time and rotate it
    quaternion_P1 = LFT1.get_quaternion(single_vector2 , vector_to_rotate_from=P_OP2)
    link = LFT1.rotate_with_quaternion(quaternion_P1, link_to_origin)

    # Move the rotated linker back to the calculated position of the phosphorus atom
    link = link + p_to_origin

    # Rotate the linekr a second time, but now the vector P_OP2 is the direction axis
    # Since we have rotated P -> OP2, we need to override the vector again from the array 'link' we just overrided
    v4 = link[id_OP2]

    # Dihedral C5' - O5' - P - OP1
    id_OP1 = LFT2.retrieve_atom_index(linker, "OP1")
    v5 = link[id_OP1]

    # Generate vector we want to rotate P_OP1 on to
    single_vector3 = generate_vector_of_interest(linker.get_OPO1(), linker.get_OP1_dihedral(), [v3, v2, v1])

    ## Rotate linker a second time
    # Generate the vectors to rotate about or to rotate with
    P_OP2 = LFT1.return_normalized(v4 - v3)  # OP2 - P ; P -> OP2
    P_OP1 = LFT1.return_normalized(v5 - v3)  # OP1 - P ; P -> OP1

    # move linker back to origin, since we do a rotation and NOT a revelation
    link = link - p_to_origin

    # we invert P_OP2, generate the quaternion and and rotate the linker
    quaternion_P2 = LFT1.get_quaternion_custom_axis(single_vector3, P_OP1, P_OP2 * -1.0)
    link = LFT1.rotate_with_quaternion(quaternion_P2, link)
    # Bring linker to the location of the phosphorus
    link = link + p_to_origin

    # Stack the arrays on top of each other
    nucleotide = np.vstack((link, nucleoside.array))

    return nucleotide

def position_next_nucleotide(next_nucleoside, prev_nucleoside, prev_linker, leading_strand : np.array) -> np.array:
    """ This function is used after position_phosphate_linker().
    It serves the purpose of adding the next nucleotide onto the leading_strand.

    next_nucleoside : json object
    leading_strand : the nucleic acid strand to which we append next_nucleoside to
    """
    # The first thing to do is to find the location of the subsequent atom, here O3', then rotate the next_nucleoside by zeta and epsilon
    # dihedral C5' - O5' - P - O3'

    shape_linker = prev_linker.get_shape()
    id_v1 = LFT2.retrieve_atom_index(prev_nucleoside, "O5'") + shape_linker
    v1 = leading_strand[id_O5]

    id_v0 = LFT2.retrieve_atom_index(prev_nucleoside, "C5'") + shape_linker
    v0 = leading_strand[id_C5]

    id_v2 = LFT2.retrieve_atom_index(prev_linker,"P")]
    v2 = leading_strand[id_v2]

    # We position the nextnuc by the position of O3'
    single_vector1 = LabF.generate_vector_of_interest(prev_linker.get_O3PO5(), next_nucleoside.get_alpha(), [v2, v1, v0])

    # Add single_vector4 to , so that we define the location of O3'
    # Multiply the normalised vector by 1.6, since 1.6 aengstrom is the distance P -> O3'
    v3 = LFT1.move_vector_to_loc(LFT1.return_normalized(single_vector1) * 1.6, v2)

    # Retrieve the vector of the newly created v3 and bring the next_nucleoside from v3_old to v3
    id_O3 = LFT2.retrieve_atom_index(next_nucleoside, "O3'")
    # Get distance from nextnuc O3' to the position defined as O3'
    O3_distance_from_origin = nextnuc.array[id_O3]
    nextnuc_distance = v6 - O3_distance_from_origin
    # Get the distance from the P_O3 now and add it to nextnuc_origin
    #   which will move the entire molecule(nextnuc) to the position of O3'
    nextnuc_loc = LFT1.move_vector_to_loc(nextnuc.array, nextnuc_distance)

    #------------ Rotate nextnuc correctly along zeta
    # O5' (v2), P (v3), O3' (v6), C3'
    P_O3_C3 = 119.032 * (np.pi / 180)
    zeta_dihr = nucleoside.get_zeta()

    v2 = nucleotide[LFT2.retrieve_atom_index(nucleoside, "O5'") + linker.get_shape()]
    single_vector5 = LabF.generate_vector_of_interest(P_O3_C3, zeta_dihr, [v6, v_P, d2])

    # We have our vector, which is O3' -> C3', so now we rotate the the nextnuc onto it
    # Get nextnuc's O3' atom to be in origin
    distance_to_origin = nextnuc_loc[id_O3]
    nextnuc_loc_tmp = LFT1.move_vector_to_origin(nextnuc_loc, distance_to_origin)

    # rotate the vector O3' -> C3' onto single_vector5
    id_C3 = LFT2.retrieve_atom_index(nextnuc, "C3'")
    v7 = nextnuc_loc[id_C3]
    O3_C3 = LFT1.return_normalized(v7 - v6)          # C3' - O3' gets the direction of O3' -> C3'
    quaternion_zeta = LFT1.get_quaternion(single_vector5, O3_C3)

    # Rotate nextnuc and move it into place
    nextnuc_loc = LFT1.rotate_with_quaternion(quaternion_zeta, nextnuc_loc_tmp)
    nextnuc_loc = nextnuc_loc + distance_to_origin

    #------------ Rotate nextnuc correctly along epsilon
    # P(v3) - O3'(v6) - C3'(v7) - C2' dihedral

    # override the current C3' vector, since we changed its location
    v7 = nextnuc_loc[id_C3]

    # Get epsilon dihedral
    epsilon_dihr = nucleoside.get_epsilon()
    # Get angle O3' - C3' - C4'
    O3_C3_C4 = 111.919 * (np.pi/180)

    single_vector6 = LabF.generate_vector_of_interest(O3_C3_C4, epsilon_dihr, [v7, v6, v_P])

    ## now that we have the vector, rotate the nucleoside appropriately.
    # get C3' and use it as the distance of the atom of interest to the origin
    C3_vector = nextnuc_loc[id_C3]
    # get C4'
    C4_vector = nextnuc_loc[LFT2.retrieve_atom_index(nextnuc, "C4'")]
    nuc_at_OG = LFT1.move_vector_to_origin(nextnuc_loc, C3_vector)
    nuc_vector_C3C4 = LFT1.return_normalized(C4_vector - C3_vector)                 # C4' - C3' gives C3' -> C4'

    O3_C3 = LFT1.return_normalized(nuc_at_OG[id_C3] - nuc_at_OG[id_O3])             # needs to be in O3' -> C3' for some reason

    # We need to rotate around the direction of the vector O3' -> C3'
    quaternion_epsilon = LFT1.get_quaternion_custom_axis(single_vector6, nuc_vector_C3C4,rotation_axis = O3_C3)
    rotated_nucatOG = LFT1.rotate_with_quaternion(quaternion_epsilon, nuc_at_OG)
    nextnuc_loc = rotated_nucatOG + C3_vector


def create_PDB_from_matrix(matrix : np.ndarray, list_of_sequence : list) -> None:
    print("Writing to pdb ...")

    df = pd.DataFrame()

    df['RecName'] = ['ATOM' for x in range(matrix.shape[0])]
    df['AtomNum'] = np.arange(start=1, stop=matrix.shape[0] + 1)
    df['AtomName'] = LFT2.pdb_AtomNames_or_ElementSymbol(list_of_sequence, "Atoms")
    df['AltLoc'] = ' '
    df['ResName'] = LFT2.pdb_Residuename(list_of_sequence)
    df['Chain'] = 'A'
    df['Sequence'] = LFT2.pdb_Sequence(list_of_sequence)
    df['X_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,0]))
    df['Y_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,1]))
    df['Z_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,2]))
    df['Occupancy'] = '1.00'
    df['Temp'] = '0.00'
    df['SegmentID'] = str('   ')
    df['ElementSymbol'] = LFT2.pdb_AtomNames_or_ElementSymbol(list_of_sequence, "Symbol")

    with open('testing_daedalus.pdb' ,'w') as pdb:
        for index, row in df.iterrows():
            split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
            pdb.write('%-6s%5s%5s%s%s%3s%5s  %8s%8s%8s%6s%6s%4s      %2s\n' % tuple(split_line))
        pdb.write('END')

