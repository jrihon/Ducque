import numpy as np
import pandas as pd
import json, sys
import time

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

        self.array =  np.asarray(json.loads(self.jason["pdb_properties"]["Coordinates"]), dtype=float)
        self.atom_list = json.loads(self.jason["pdb_properties"]["Atoms"])

    # One method that makes it easy to call the prompted dihedrals of the nucleoside.
    # because the dihedral is still inside a the backbone dictionary, we need to load the string (json.loads) again
    def get_dihedral(self, dihedral : str) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])[dihedral])

    def get_angle(self, angle : str) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])[angle]) * (np.pi/180)


class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_shape(self) -> int:
        # returns the size of the linker shape, but only the first value
        return int(json.loads(self.jason["pdb_properties"]["Shape"])[0])


    def get_COP(self) -> float:
        # since the angles are inside the angles dictionary, we don"t need to load string again
        return float(self.jason["Angles"]["C5_O5_P"]) * (np.pi / 180)

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
    single_vector = LFT1.generate_and_rotate_single_vector(theta_interpolate, angle, quaternion)

    return single_vector

def position_phosphate_linker(nucleoside, nucl_array : np.ndarray, linker) -> np.ndarray:
    """
    This function positions and rotates the phosphate linker onto the nucleoside of interest.
    First the location of the phosphorus atom is calculated for.
    Secondly, the first rotation is carried out based on the P -> OP2 vector.
    Lastly, the second rotation is carried out based on the P -> OP1 vector.

    nucleoside is a json object
    linker is a json object
    """
    # Atom parsing list : last ones of the nucleoside and the ones needed from the linker
    ATP = ["C4'", "C5'", "O5'", "P", "OP2", "OP1"]

    # Retrieve the vectors of the atoms that make up the dihedral you research
    # Dihedral C4' - C5' - O5' - P
    id_v0 = LFT2.retrieve_atom_index(nucleoside, ATP[0])
    v0 = nucl_array[id_v0]

    id_v1 = LFT2.retrieve_atom_index(nucleoside, ATP[1])
    v1 = nucl_array[id_v1]

    id_v2 = LFT2.retrieve_atom_index(nucleoside, ATP[2])
    v2 = nucl_array[id_v2]

    # Find the vector that corresponds to the O5' -> P vector
    single_vector1 = generate_vector_of_interest(linker.get_COP(), nucleoside.get_dihedral("beta"), [v2, v1, v0])

    # Add the single_vector1 to O5' atom, to denote the location of P. Then move the linker array to that location
    link = LFT1.position_phosphate(v2, single_vector1, linker.array)   # the linker is in the correct position  

    ## We will rotate the linker twice, to get the correct orientation of the linker in 3D space
    # Dihedral C5' - O5' - P - OP2
    id_v3 = LFT2.retrieve_atom_index(linker, ATP[3])
    v3 = link[id_v3]

    single_vector2 = generate_vector_of_interest(linker.get_OPO2(), linker.get_OP2_dihedral(), [v3, v2, v1])

    # This is the distance from the linker's atom to the origin
    link_distance = link[id_v3]

    # move the linker to the origin, by positioning the phosphorus at [0,0,0]
    link_to_origin = LFT1.move_vector_to_origin(link, link_distance)

    # Rotate the linker a first time
    # Define the vector that goes from P to OP2 and normalize it
    id_v4 = LFT2.retrieve_atom_index(linker, ATP[4])
    v4 = link[id_v4]
    p3_4 = LFT1.return_normalized(v4 - v3)

    # Get quaternion to rotate the linker a first time and rotate it
    quaternion_P1 = LFT1.get_quaternion(single_vector2 , p3_4)
    link = LFT1.rotate_with_quaternion(quaternion_P1, link_to_origin)

    # Move the rotated linker back to the calculated position of the phosphorus atom
    link = LFT1.move_vector_to_loc(link, link_distance)

    # Rotate the linekr a second time, but now the vector P_OP2 is the direction axis
    # Since we have rotated the linker, we need to override the vector again from the array 'link' we just overrided
    v4 = link[id_v4]

    # Dihedral C5' - O5' - P - OP1
    id_v5 = LFT2.retrieve_atom_index(linker, ATP[5])
    v5 = link[id_v5]

    # Generate vector we want to rotate P_OP1 on to
    single_vector3 = generate_vector_of_interest(linker.get_OPO1(), linker.get_OP1_dihedral(), [v3, v2, v1])

    # normalize the single vector, multiply with the set distance (P-O) and replace it with the index of OP1 in the link array, making it the new vector for OP1
    distance_P_O = 1.48
    link[id_v5] = (LFT1.return_normalized(single_vector3) * distance_P_O) + link[id_v3]

    # Stack the arrays on top of each other
    nucleotide = np.vstack((link, nucl_array))

    return nucleotide

def position_next_nucleoside(next_nucleoside, prev_nucleoside, prev_linker, leading_strand : np.array) -> np.array:
    """ This function is used after position_phosphate_linker().
    It serves the purpose of adding the next nucleotide onto the leading_strand.

    next_nucleoside : json object
    leading_strand : the nucleic acid strand to which we append next_nucleoside to
    """
    # The first thing to do is to find the location of the subsequent atom, here O3', then rotate the next_nucleoside by zeta and epsilon
    # Afterwards, we turn over the epsilon dihedral and by rotate the normal of the plane have to the plane we want; this positions everything!

    # Atom Parsing List (ATP) = Parse which linker and which nucleotide the previous one is
    APL = ["C5'", "O5'", "P", "O3'", "C3'", "C4'"]

    # Dihedral Parsing List (DPL) = Parse which dihedrals are required to rotate on and over
    DPL = ["alpha", "zeta", "epsilon"]

    # Angle Parsing List (AngPL)
    AngPL = [prev_linker.get_O3PO5(), 119.032 * (np.pi / 180), ]

    # List of which we parse our 
    #### POSITION THE NEXT NUCLEOSIDE PROPERLY TO O3'
    # dihedral C5' - O5' - P - O3'
    # get the size of the linker, since the array goes : " linker.array + nucleotide.array", where linker.array precedes the indices of the nucleotide.array
    shape_linker = prev_linker.get_shape()
    id_v0 = LFT2.retrieve_atom_index(prev_nucleoside, APL[0]) + shape_linker
    v0 = leading_strand[id_v0]

    id_v1 = LFT2.retrieve_atom_index(prev_nucleoside, APL[1]) + shape_linker
    v1 = leading_strand[id_v1]

    id_v2 = LFT2.retrieve_atom_index(prev_linker,APL[2])
    v2 = leading_strand[id_v2]

    # Get angle and dihedral:
    alpha_angle = AngPL[0]
    alpha_dihr = next_nucleoside.get_dihedral(DPL[0])

    ## We position the nextnuc by the position of O3'
    single_vector1 = generate_vector_of_interest(alpha_angle, alpha_dihr, [v2, v1, v0])
    linker_nuc_distance = 1.6
    # Add single_vector1 to v2, so that we define the location of ATP[3]
    v3 = LFT1.move_vector_to_loc(LFT1.return_normalized(single_vector1) * linker_nuc_distance, v2)

    ## Bring the nucleoside array to the new position
    id_v3 = LFT2.retrieve_atom_index(next_nucleoside, APL[3])
    v3_original = next_nucleoside.array[id_v3]
    # Get distance from next_nucleoside O3' (v3_old) to the position defined as O3' (v3), so the distance that we need to move the array with
    distance_v3_v3_original = v3 - v3_original
    # Get the distance from the P -> O3' now and add it to nextnuc_origin which will move the entire molecule(nextnuc) to the position of O3'
    next_nucleoside_loc = LFT1.move_vector_to_loc(next_nucleoside.array, distance_v3_v3_original)

#    for i in range(1, len(DTP)):
#        #Do the final vector rotation + plane rotation
#        # Only at the final part do we need to do the rotation.
#        # We only rotate over the previous dihedral and then the following dihedral.
#        # The reason we forloop is because perhaps there might be cases where just a single rotation and a plane rotation are sufficient, like with morpholino's presumably.
#        if i + 1 == len(DTP):
#            pass
#            # return
#
#        # Just continue if it is not the last one in the loop
#        id_vA = LFT2.retrieve_atom_index(prev_nucleoside, APL[i]) + shape_linker
#        vA = leading_strand[id_vA]
#
#        id_vB = LFT2.retrieve_atom_index(linker, APL[i + 1])
#        vB = leading_strand[id_vB]
#
#        id_vC = LFT2.retrieve_atom_index(next_nucleoside, APL[i + 2])
#        vC = next_nucleoside_loc[id_vC]
#
#        dihedral_N = next_nucleoside.get_dihedral(DPL[i])
#        angle_N = AngPL[i]
#
#        single_vector_N = generate_cone_vector(angle_N, dihedral_N, [vC, vB, vA]
    #### POSITION THE NUCLEOSIDE ALONG THE EPSILON BOND BY ROTATIONG THE ZETA DIHEDRAL
    ## Now, we define the next dihedral and find the atom. This is likely to be zeta, so we rotate over zeta afterwards
    # O5' (v1), P (v2), O3' (v3), C3'
    zeta_angle = 119.032 * (np.pi / 180)
    zeta_dihr = next_nucleoside.get_dihedral("zeta")

    single_vector2 = generate_vector_of_interest(zeta_angle, zeta_dihr, [v3, v2, v1])

    # We have our vector so now we rotate next_nucleoside onto single_vector2
    id_v4 = LFT2.retrieve_atom_index(next_nucleoside, ATP[4])
    v4 = next_nucleoside_loc[id_v4]
    p3_4 = LFT1.return_normalized(v4 - v3)
    quaternion_zeta = LFT1.get_quaternion(single_vector2, p3_4)

    ## Move next_nucleoside to the origin, rotate it and move it back into place
    # Get nextnuc's O3' atom to be in origin
    distance_to_origin = next_nucleoside_loc[id_v3]

    next_nucleoside_originloc = LFT1.move_vector_to_origin(next_nucleoside_loc, distance_to_origin)
    next_nucleoside_originloc = LFT1.rotate_with_quaternion(quaternion_zeta, next_nucleoside_originloc)
    next_nucleoside_loc = LFT1.move_vector_to_loc(next_nucleoside_originloc, distance_to_origin)

    ## NOW ROTATE THE PLANE ALONG THE NEXT DIHEDRAL, BUT DO NOT OVERRIDE THE ONES AT ARE ALREADY IN PLACE
    # First we rotate the nucleoside as usual, now on epsilon
    epsilon_dihr = next_nucleoside.get_dihedral("epsilon")
    epsilon_angle = 111.919 * (np.pi/180)
    # override vectors that have moved (technically v3 has not moved but whatever, just in case)
    v4 = next_nucleoside_loc[id_v4] ; v3 = next_nucleoside_loc[id_v3]
    # retrieve last vector.
    id_v5 = LFT2.retrieve_atom_index(next_nucleoside, ATP[5])
    v5 = next_nucleoside_loc[id_v5]
    # Search for the single_vector
    single_vector_n = generate_vector_of_interest(epsilon_angle, epsilon_dihr, [v4, v3, v2])

    # move the nucleoside to the origin
    distance_to_origin_n = next_nucleoside_loc[id_v4]
    next_nucleoside_originloc = LFT1.move_vector_to_origin(next_nucleoside_loc, distance_to_origin_n)

    ## Generate to normal vectors of the planes of interest. The order in which you perform the cross product is not important BUT!!!
    # It IS important that the same vectors are operated on in the same order for both normal vectors!!!
    # normal to rotate from
    n0 = LFT1.get_normal_vector_of_plane(v4 - v3, v5 - v4)
    # normal to rotate to
    n1 = LFT1.get_normal_vector_of_plane(v4 - v3, single_vector_n)

    # Get the indices of the vectors you do not want moved
    exclusion_list = ["O3'", "C3'"]
    exclude_these_atoms = [LFT2.retrieve_atom_index(next_nucleoside, atom) for atom in exclusion_list]

    # The order of the quaternion does matter, as it starts with " vector to rotate to" and secondly with "vector that we want to rotate from "
    quaternion_plane = LFT1.get_quaternion(n1, n0)
    next_nucleoside_originloc_rotated = LFT1.apply_rotation_of_planes(quaternion_plane, next_nucleoside_originloc, exclude_these_atoms)
    next_nucleoside_loc = next_nucleoside_originloc_rotated + distance_to_origin_n

    return next_nucleoside_loc


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

    filename = "testing_daedalus.pdb"
    with open(filename ,'w') as pdb:
        for index, row in df.iterrows():
            split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
            pdb.write('%-6s%5s%5s%s%s%3s%5s  %8s%8s%8s%6s%6s%4s      %2s\n' % tuple(split_line))
        pdb.write('END')

