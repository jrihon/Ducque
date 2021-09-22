import numpy as np
import pandas as pd
import json, sys
from typing import Union

import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2
import labyrinth_func_tools3 as LFT3
import labyrinth_func_andthefunkybunch as LFFB


""" Labyrinth_func.py
The script that contains the classes and all the functions that concatenate the workflow of consecutively adding the linker and nucleotides. """

                                                                             #### CLASSES
class Nucleoside:
    """ nucleoside = Nucleoside(codex_acidum_nucleicum[nucleic_acid-string][0]) """

    def __init__(self, jsonfile):

        with open(jsonfile, "r") as jsonf:
            self.jason = json.load(jsonf)

        self.array =  np.asarray(json.loads(self.jason["pdb_properties"]["Coordinates"]), dtype=float)
        self.atom_list = json.loads(self.jason["pdb_properties"]["Atoms"])
        self.mol_length = int(json.loads(self.jason["pdb_properties"]["Shape"])[0])

    def get_dihedral(self, dihedral : str) -> float:                        # because the dihedral is still inside a dictionary, we need to load the string (json.loads)
        """ return dihedral value of the queried dihedral """
        return float(json.loads(self.jason["angles"]["dihedrals"])[dihedral])

    def get_angle(self, angle : str) -> float:
        """ return angle value of the queried angle. Needs to be converted to radians """
        return float(json.loads(self.jason["angles"]["bond_angles"])[angle]) * (np.pi/180)

    def get_base_denominator(self) -> str:
        """ returns the type of base of the nucleic acid. So if the base is Guanosine, return 'G'. """
        return json.loads(self.jason["identity"])[2][-1]

    def get_nucleic_acid_code(self) -> str:
        """ returns the type of chemistry of the nucleic acid. """
        return json.loads(self.jason["identity"])[1]


class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_shape(self) -> int:
        """ returns the size of the linker shape, but only the first value """
        return int(json.loads(self.jason["pdb_properties"]["Shape"])[0])

    def get_nucleic_acid_code(self) -> str:
        """ returns the type of chemistry of the nucleic acid. """
        return json.loads(self.jason["identity"])[0]



                                                                             #### FUNCTIONS
# Functions that are meant to concatenate and bypass the iterative coding in Labyrinth.py : Architecture() 
def generate_vector_of_interest(angle : float, dihedral : float, atom_array : np.array) -> np.array:
    """ This function generates a single vector for which there exists only one angle and dihedral.
    The atom_array contains the three first atoms in the sequence that make up the dihedral.
    Example = if the sequence of a dihedral is C4' - C5' - O5' - P (beta backbone), then the atom_array is [O5', C5', C4'] """

    # Get the vector to rotate the cone vector onto. This is done on the middle two atoms of the sequence
    # Example : [O5'] minus (-) [C5'] results in a vector that goes [C5' -> O5']
    p0 = LFT1.return_normalized(atom_array[0] - atom_array[1])

    # Generate the cone vector
    cone_vector = LFT1.generate_cone_vector(angle)

    # Get the quaternion that corresponds to the desired rotation
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
    APL = LFT3.Atom_Parsing_List(nucleoside, linker)

    # Retrieve the vectors of the atoms that make up the dihedral you research
    # Dihedral C4' - C5' - O5' - P
    id_v0 = LFT2.retrieve_atom_index(nucleoside, APL[0])
    id_v1 = LFT2.retrieve_atom_index(nucleoside, APL[1])
    id_v2 = LFT2.retrieve_atom_index(nucleoside, APL[2])
    id_v3 = LFT2.retrieve_atom_index(linker, APL[3])

    v0 = nucl_array[id_v0]
    v1 = nucl_array[id_v1]
    v2 = nucl_array[id_v2]

    # Find the vector that corresponds to the O5' -> P vector
    single_vector1 = generate_vector_of_interest(nucleoside.get_angle("beta"), nucleoside.get_dihedral("beta"), [v2, v1, v0])

    # Position the phosphate by moving the linker's center atom to the origin and then bringing to the location that was calculated for
    move_link_to = v2 + (single_vector1 * 1.6)
    link = LFT1.move_vector_to_loc(linker.array - linker.array[id_v3], move_link_to)

    ## We will rotate the linker twice, to get the correct orientation of the linker in 3D space
    # Dihedral C5' - O5' - P - OP2
    v3 = link[id_v3]

    single_vector2 = generate_vector_of_interest(linker.get_angle("OPO"), linker.get_dihedral("OP2_dihedral"), [v3, v2, v1])

    # This is the distance from the linker's atom to the origin
    link_distance = link[id_v3]

    # move the linker to the origin, by positioning the phosphorus at [0,0,0]
    link_to_origin = LFT1.move_vector_to_origin(link, link_distance)

    # Rotate the linker a first time
    # Define the vector that goes from P to OP2 and normalize it
    id_v4 = LFT2.retrieve_atom_index(linker, APL[4])
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
    id_v5 = LFT2.retrieve_atom_index(linker, APL[5])
    v5 = link[id_v5]

    # Generate vector we want to rotate P_OP1 on to
    single_vector3 = generate_vector_of_interest(linker.get_angle("OPO"), linker.get_dihedral("OP1_dihedral"), [v3, v2, v1])

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
        leading_strand : the nucleic acid strand to which we append next_nucleoside to. """
    # The first thing to do is to find the location of the subsequent atom, here O3', then rotate the next_nucleoside by zeta and epsilon
    # Afterwards, we turn over the epsilon dihedral and by rotate the normal of the plane have to the plane we want; this positions everything!

    # Atom Parsing List (ATP) = Parse which linker and which nucleotide the previous one is
    APL = LFT3.Atom_Parsing_List(prev_nucleoside, prev_linker, next_nucleoside)
    # Dihedral Parsing List (DPL) = Parse which dihedrals are required to rotate on and over
    DPL = ["alpha", "zeta", "epsilon"]
    # Angle Parsing List (AngPL)
    AngPL = ["alpha", "zeta", "epsilon"]

    #### POSITION THE NEXT NUCLEOSIDE PROPERLY.
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
    alpha_angle = next_nucleoside.get_angle(AngPL[0])
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

    for i in range(1, len(DPL)):
        # Only at the final dihedral do we require a plane rotation.
        # The reason we forloop is because perhaps there might be cases where just a single rotation and a plane rotation are sufficient, like with morpholino's presumably.
        if i + 1 == len(DPL):
            # Do the final vector rotation + plane rotation
            id_vA = LFT2.retrieve_atom_index(prev_linker, APL[i])
            vA = leading_strand[id_vA]
            id_vB = LFT2.retrieve_atom_index(next_nucleoside, APL[i + 1])
            vB = next_nucleoside_loc[id_vB]
            id_vC = LFT2.retrieve_atom_index(next_nucleoside, APL[i + 2])
            vC = next_nucleoside_loc[id_vC]
            #Get the required angles
            angle_N = next_nucleoside.get_angle(AngPL[i])
            dihedral_N = next_nucleoside.get_dihedral(DPL[i])

            single_vector_N = generate_vector_of_interest(angle_N, dihedral_N, [vC, vB, vA])

            # Retrieve the next atom vector in the sequence, required for rotation
            id_vD = LFT2.retrieve_atom_index(next_nucleoside, APL[i + 3])
            vD = next_nucleoside_loc[id_vD]

            ## Generate to normal vectors of the planes of interest. The order in which you perform the cross product is not important BUT!!!
            # It IS important that the same vectors are operated on in the same order for both normal vectors!!!
            # normal to rotate from
            n0 = LFT1.get_normal_vector_of_plane(vC - vB, vD - vC)
            # normal to rotate to
            n1 = LFT1.get_normal_vector_of_plane(vC - vB, single_vector_N)
            # The order of the quaternion does matter, as it starts with " vector to rotate to" and secondly with "vector that we want to rotate from "
            quaternion_plane = LFT1.get_quaternion(n1, n0)

            #Now bring the array to the origin
            distance_to_origin_N = next_nucleoside_loc[id_vC]
            next_nucleoside_loc = LFT1.move_to_origin_ROTATE_move_back_to_loc(quaternion_plane, next_nucleoside_loc, distance_to_origin_N)

            return next_nucleoside_loc

        # Just continue if it is not the last one in the loop
        id_vA = LFT2.retrieve_atom_index(prev_nucleoside, APL[i]) + shape_linker
        vA = leading_strand[id_vA]
        id_vB = LFT2.retrieve_atom_index(prev_linker, APL[i + 1])
        vB = leading_strand[id_vB]
        id_vC = LFT2.retrieve_atom_index(next_nucleoside, APL[i + 2])
        vC = next_nucleoside_loc[id_vC]

        #Get the required angles
        angle_N = next_nucleoside.get_angle(AngPL[i])
        dihedral_N = next_nucleoside.get_dihedral(DPL[i])

        single_vector_N = generate_vector_of_interest(angle_N, dihedral_N, [vC, vB, vA])

        # Retrieve the appropriate quaternion for the rotation 
        id_vD = LFT2.retrieve_atom_index(next_nucleoside, APL[i + 3])
        vD = next_nucleoside_loc[id_vD]
        pC_D = LFT1.return_normalized(vD - vC)
        quaternion_N = LFT1.get_quaternion(single_vector_N, pC_D)

        #Now rotate your molecule onto single_vector_N
        distance_to_origin_N = next_nucleoside_loc[id_vC]

        next_nucleoside_loc = LFT1.move_to_origin_ROTATE_move_back_to_loc(quaternion_N, next_nucleoside_loc, distance_to_origin_N)


def assert_leading_strand_nucleotide_conformation(conformations, prev_nucleoside, prev_link, leading_strand):
    """ Position the different conformations and assert which conformation is the most parralel with respect to the bases of both nucleosides """ 

    # If there is only one conformation possible, no point in asserting the best one
    if len(conformations) == 1:
        return conformations[0]


    # Position the different conformations
    arrays_of_the_conformations = np.zeros(len(conformations), dtype=object)

    for conf in range(len(conformations)):
        nucleoside = Nucleoside(conformations[conf])
        arrays_of_the_conformations[conf] = position_next_nucleoside(nucleoside, prev_nucleoside, prev_link, leading_strand)

    # Get the json object parse the correct atoms for the vector comparison
    next_nucleoside = Nucleoside(conformations[0])

    prev_base = prev_nucleoside.get_base_denominator() 
    next_base = next_nucleoside.get_base_denominator()
    prev_base_atoms, next_base_atoms = LFT3.retrieve_atoms_for_plane_rotation_of_complement(prev_base, next_base) 

    # Parse the indexes of the atoms of interest
    prev_base_indexes = LFT2.retrieve_atom_index_MULTIPLE(prev_nucleoside, prev_base_atoms, index_counter=prev_link.mol_length)
    next_base_indexes = LFT2.retrieve_atom_index_MULTIPLE(next_nucleoside, next_base_atoms)

    # Previous nucleoside Cross Vector
    cross_prev = LFFB.return_cross_vector_for_plane_rotation(leading_strand, prev_base_indexes)

    # Array of different cross vectors of the next nucleoside
    possible_cross_vectors = np.zeros(len(conformations), dtype=object)

    for conf in range(len(conformations)):
        # Previous nucleoside Cross Vector
        nucleoside = Nucleoside(conformations[conf])
        possible_cross_vectors[conf] = LFFB.return_cross_vector_for_plane_rotation(arrays_of_the_conformations[conf], next_base_indexes)

    # Compare the dot product of the vector products of the next nuc with that of the previous nuc
    scalar_products_array = np.zeros(len(conformations), dtype=object)

    for v in range(len(conformations)):
        scalar_products_array = np.dot(cross_prev, possible_cross_vectors[v])

    smallest_val = scalar_products_array.min()

    index_smallest_val = np.where(scalar_products_array == smallest_val)[0]

    return conformations[index_smallest_val]


def position_complementary_base(leading_base, complementary_base, leading_array : np.ndarray, index_lead : int) -> np.ndarray:
    """ This functions positions the base correctly, which inherently positions the entire nucleoside correctly.
        After this function, there comes an iterative process of fitting the backbone correctly.

        leading_base is a json object
        complementary_base is a json object                     """

    ### DATA PARSING
    ## Parse the correct atoms to get started for
    # Parse the type of base from the nucleoside json object
    leadingBase = leading_base.get_base_denominator()
    complBase = complementary_base.get_base_denominator()
    rotPlane_leadingBase_atoms, rotPlane_complBase_atoms = LFT3.retrieve_atoms_for_plane_rotation_of_complement(leadingBase, complBase)
    ## Get the proper vectors for the plane rotation
    # Leading base vector
    rotPlane_leadingBase_atoms_id = LFT2.retrieve_atom_index_MULTIPLE(leading_base, rotPlane_leadingBase_atoms, index_lead)
    vCross_leadingBase = LFFB.return_cross_vector_for_plane_rotation(leading_array, rotPlane_leadingBase_atoms_id)

    # Complementary base vector
    rotPlane_complBase_atoms_id = LFT2.retrieve_atom_index_MULTIPLE(complementary_base, rotPlane_complBase_atoms)
    vCross_complBase = LFFB.return_cross_vector_for_plane_rotation(complementary_base.array, rotPlane_complBase_atoms_id)
    
    ### ROTATE THE PLANE
    ## Rotate the complementary base onto the plane of the leading base
    # Get the cross product of the atoms of leading_base that make it if both bases are in the same plane, the cross product vectors are exactly opposite
    plane1_quaternion = LFT1.get_quaternion(vCross_leadingBase * -1.0, vCross_complBase)

    ## Apply the rotation
    # vC0 is atom that binds the base with the sugar moiety, making it the central atom over which we rotate when we bring the molecule to the origin
    vC0 = complementary_base.array[rotPlane_complBase_atoms_id[0]]
    tmp_compl_arr = LFT1.move_vector_to_origin(complementary_base.array, vC0)
    compl_nucleoside_array = LFT1.rotate_with_quaternion(plane1_quaternion, tmp_compl_arr)

    # Get the required parameters for the rotations
    Q_angle, Q_dihedral, R_angle, R_dihedral = LFT3.retrieve_angles_and_dihedrals_for_initial_base_positioning(leadingBase)
    Q_dist_between_bases = 1.81
    R_dist_between_bases = 1.87

    ## Apply a translation to get the initial positioning
    Q_leadingBase_atoms, Q_complBase_atom, Q_distance = LFT3.retrieve_atoms_for_positioning_of_complement1(leadingBase, complBase)
    Q_leadingBase_atoms_id = LFT2.retrieve_atom_index_MULTIPLE(leading_base, Q_leadingBase_atoms, index_lead)

    # Calculate the position at which we want to the R positioning to be. See paper @Figure X to show what R positioning is.
    Q_position = LFFB.return_position_for_complementary_base(Q_angle, Q_dihedral, leading_array, Q_leadingBase_atoms_id, Q_distance)

    ## Now that we have the Q_position, we move the complementary base to the position by doing the following :
    # Get distance at which we want to position the base's hydrogen bond
    Q_vector = leading_array[Q_leadingBase_atoms_id[2]] + Q_position
    # retrieve the vector of the atom that we want to attach 'Q_vector' to.
    Q_complBase_atom_id = LFT2.retrieve_atom_index(complementary_base, Q_complBase_atom)

    # the distance from the location we want to have compl1_base2_atom move to. But we move the array from where it is at now, so substract that and then move it
    move_to_Q = Q_vector - compl_nucleoside_array[Q_complBase_atom_id]

    compl_nucleoside_array = LFT1.move_vector_to_loc(compl_nucleoside_array, move_to_Q)

    ### After positioning the planes and moving the nucleoside to Q, we now rotate the nucleoside's base to R
    R_leadingBase_atoms, R_complBase_atom, R_distance = LFT3.retrieve_atoms_for_position_of_complement2(leadingBase, complBase)
    R_leadingBase_atoms_id = LFT2.retrieve_atom_index_MULTIPLE(leading_base, R_leadingBase_atoms, index_lead)

    R_position = LFFB.return_position_for_complementary_base(R_angle, R_dihedral, leading_array, R_leadingBase_atoms_id, R_distance)

    # Creates the coordinate at which we want to move our the Q_complement atom, from the complementary nucleoside's base to
    R_move_to = leading_array[R_leadingBase_atoms_id[2]] + R_position
    R_compl_id = LFT2.retrieve_atom_index(complementary_base, R_complBase_atom)
    R_vector = compl_nucleoside_array[Q_complBase_atom_id]

    # Create quaternion to rotate the complementary base a final time
    j1 = LFT1.return_normalized(R_move_to - R_vector)
    j2 = LFT1.return_normalized(compl_nucleoside_array[R_compl_id] - R_vector)

    plane2_quaternion = LFT1.get_quaternion(j1, j2)

    return LFT1.move_to_origin_ROTATE_move_back_to_loc(plane2_quaternion, compl_nucleoside_array, R_vector)


def tilt_array_to_get_a_better_fit(compl_nuc, compl_linker, prev_nuc, prev_linker, compl_nuc_arr : np.ndarray, complementary_strand : np.ndarray, index_compl : int) -> np.ndarray :
    """ Have the returned array from assert_possible_base_conformations fit the leading strand better. Here we will turn the return array according to length, then dihedral
        compl_nuc is a json object
        prev_nuc is a json object   """

    ## Now we tilt the array slighted if necessary
    # The 'if necessary' will be depending on the length and the dihedral it makes

    # Atom Parsing List for the knowing which atoms to parse from the respective arrays ; for bond length and dihedral evaluation
    APL = LFT3.Atom_Parsing_List(compl_nuc, compl_linker, prev_nuc)
    dihedral = compl_nuc.get_dihedral("alpha")
    angle = compl_nuc.get_angle("alpha")

    # Parse the last atom needed from the complementary strand. NB : index_compl has been set to the correct value in the function ' assert_possible ... _and_fit() '.
    # This is the last final atom in the backbone ( 3'-> 5), so in the case of DNA, this would store the vector of the phosphorus atom in its linker
    id_compl_strand = LFT2.retrieve_atom_index(prev_nuc, APL[3]) + index_compl + prev_linker.mol_length
    v_compl_strand = complementary_strand[id_compl_strand]

    # Parse the other atoms needed from the current nucleotide that is being fitted
    id_v0 = LFT2.retrieve_atom_index(compl_nuc, APL[0]) + compl_linker.mol_length
    id_v1 = LFT2.retrieve_atom_index(compl_nuc, APL[1]) + compl_linker.mol_length
    id_v2 = LFT2.retrieve_atom_index(compl_linker, APL[2])

    v0 = compl_nuc_arr[id_v0]
    v1 = compl_nuc_arr[id_v1]
    v2 = compl_nuc_arr[id_v2]

    # Calculate of the distance. We will check later if it is between 1 and 2 Angstrom.
    calculated_length = LFT1.get_length_of_vector(v2, v_compl_strand)

    # Calculate the dihedral. Let's say the dihedral should not deviate more than 25 degrees?
    calculated_dihr = LFT1.dihedral_single(v0, v1, v2, v_compl_strand)

    # Create booleans to assert the calculated values and see if they are within ranges of concordance
    bool_length = LFT1.assert_length_of_vector(calculated_length)

    bool_dihedral = LFT1.assert_dihedral_of_nucleotide(calculated_dihr, dihedral)

    # If both values are true, just return the array as is, since it is already relatively well positioned
    if bool_length and bool_dihedral :
        return compl_nuc_arr

    # If one or both values are false, try to position the bases to an appropriate orientation
    else :
        ## Generate a quaternion that orients the alpha dihedral properly
        angle_of_rot_TEST = 5 * (np.pi/180)
        # Retrieve which the denominator of the base that the nucleotide is
        compl_base = compl_nuc.get_base_denominator()

        # Make a custom direction and angle. Not to worry, the direction is perpendicular to the movement
        compl_base_atoms1, _ = LFT3.retrieve_atoms_for_plane_rotation_of_complement(compl_base, compl_base)
        compl_base_atoms1_0 = LFT2.retrieve_atom_index(compl_nuc, compl_base_atoms1[0]) + compl_linker.mol_length

        compl_base_atoms2, _ , __= LFT3.retrieve_atoms_for_positioning_of_complement1(compl_base, compl_base)
        compl_base_atoms2_2 = LFT2.retrieve_atom_index(compl_nuc, compl_base_atoms2[2]) + compl_linker.mol_length

        v_direction = LFT1.return_normalized(compl_nuc_arr[compl_base_atoms2_2] - compl_nuc_arr[compl_base_atoms1_0])

        # Rotate both directions the nuc1 and assert the distance between nuc1 and nuc2, then remember the direction of the smallest turn
        nuc1_quaternion = LFT1.get_custom_quaternion(v_direction, angle_of_rot_TEST)
        test_turn_arr = LFT1.move_to_origin_ROTATE_move_back_to_loc(nuc1_quaternion, compl_nuc_arr, compl_nuc_arr[compl_base_atoms1_0])
        # Distance between phosphorus of the previous nucleotide and the (i.e.) O3' of the next nucleoside (in the example given a DNA nucleoside)
        dist_between_nucs = LFT1.get_length_of_vector(test_turn_arr[id_v2], complementary_strand[id_compl_strand])

        # Check if the distance between subsequent nucleotides is suitable, according to the set boundaries.
        # If that is not the case, store the current distance, rotate the nucleotide in the reverse direction and store that new distance between nucleotides too.
        # Check which distance is the smallest and remember the orientation of the direction axis, as this will be the direction to rotate over from now on.
        if not LFT1.assert_length_of_vector(dist_between_nucs):
            x1 = dist_between_nucs

            v_direction_test = v_direction * -1.0
            nuc1_quaternion = LFT1.get_custom_quaternion(v_direction_test, angle_of_rot_TEST)
            test_turn_arr = LFT1.move_to_origin_ROTATE_move_back_to_loc(nuc1_quaternion, compl_nuc_arr, compl_nuc_arr[compl_base_atoms1_0])
            x2 = LFT1.get_length_of_vector(test_turn_arr[id_v2], complementary_strand[id_compl_strand])

            if x1 > x2 :
                v_direction *= -1.0

        # Generate angles of rotation. Make note that we won't go past a 15 degree angle, since that general will not be necessary to turn that much.
        array_of_rot_angles = np.linspace(2.5, 30, 8) * (np.pi/180)

        stored_distances = np.zeros(len(array_of_rot_angles))
        for i in range(len(array_of_rot_angles)) :
            nuc1_quaternion = LFT1.get_custom_quaternion(v_direction, array_of_rot_angles[i])
            testnuc1 = LFT1.move_to_origin_ROTATE_move_back_to_loc(nuc1_quaternion, compl_nuc_arr, compl_nuc_arr[compl_base_atoms1_0])
            dist_between_nucs = LFT1.get_length_of_vector(testnuc1[id_v2], complementary_strand[id_compl_strand])
            stored_distances[i] = dist_between_nucs

            # If the distance is suitable according to the boundaries, return the nucleoside's array
            if LFT1.assert_length_of_vector(dist_between_nucs) :
                return testnuc1

        # If this part of the code is reached, that means none of the evaluated distances are suitable. Let's search the shortest one in the list
        min_val = stored_distances.min()
        index_min = np.where(stored_distances == min_val)
        angle_of_rot = array_of_rot_angles[index_min[0]]
        nuc1_quaternion = LFT1.get_custom_quaternion(v_direction, float(angle_of_rot))

        return LFT1.move_to_origin_ROTATE_move_back_to_loc(nuc1_quaternion, compl_nuc_arr, compl_nuc_arr[compl_base_atoms1_0])


def assert_starting_bases_of_complementary_strand(compl1_base_confs : list, compl2_base_confs : list, compl2_linker, lead_bases : list, leading_strand : np.ndarray, index_lead : int) -> np.ndarray :
    """ The first two bases need to be fitted as best as it can.
        Unfortunately, the first nucleoside incorporated has no reference point. As such, the second base ( a nucleotide ) cannot be fitted correctly.
        compl2_linker is a json object

        The thoughtprocess here is that we calculate the magnitude of the tilt for the second base and rotate both bases by half of that magnitude.
        So we first assert all the possibilities by iteratively fitting the bases.

        Propose : the amount of conformations is equal to 2 for both nucleosides. This sets the total amount of combinations to 4.

        conf1   |   conf1
        conf1   |   conf2
        --------|--------
        conf1   |   conf2
        conf2   |   conf2

        These combinations are asserted and the best one is selected for fitting. """

    # Instantiate an array that holds all the possible conformation combinations
    # Consider we have 3 conformations for the 1st base and 2 for the 2nd. with this we create an array that has 2 rows and 3 columns, making it easier to loop over the
    #   array based on position of the nucleotide
    array_of_possible_dinucleotide_conformations = np.zeros((len(compl1_base_confs), len(compl2_base_confs)), dtype=object)

    # Create instanced objects of the two nucleotides to parse the correct vectors
    lead1_base, lead1_link = Nucleoside(lead_bases[0][0]), Desmos(lead_bases[0][1])
    lead2_base, lead2_link = Nucleoside(lead_bases[1][0]), Desmos(lead_bases[1][1])

    # Prepare variables for the forloop
    index_lead -= lead1_base.mol_length
    size_of_lead_base2 = lead2_base.mol_length + lead2_link.mol_length

    # This works for for when there is a single conformation or multiple conformations
    for i in range(len(compl1_base_confs)):
        # Instantiate the complementary base
        compl1_base = Nucleoside(compl1_base_confs[i])
        # position the complementary base
        compl1_base_arr = position_complementary_base(lead1_base, compl1_base, leading_strand, index_lead)

        for j in range(len(compl2_base_confs)):
            #decrement the index_lead
            index_lead -= size_of_lead_base2

            # Instantiate the complementary base
            compl2_base = Nucleoside(compl2_base_confs[j])
            # position the complementary base
            compl2_base_arr = position_complementary_base(lead2_base, compl2_base, leading_strand, index_lead)
            # add linker to the complementary base
            compl2_base_arr = position_phosphate_linker(compl2_base, compl2_base_arr, compl2_linker)

            # Store the array of the dinucleotide
            array_of_possible_dinucleotide_conformations[i][j] = np.vstack((compl1_base_arr, compl2_base_arr))

            #increment the index_lead back
            index_lead += size_of_lead_base2


    # decrement it again, since we return this variable for the upcoming nucleotides after this function finalises.
    index_lead -= size_of_lead_base2

    # Parse the index of the values we want. compl1 is the atom that is being attached to the linker of the previous nucleotide. compl2 is the last atom in the backbone of the linker.
    APL1 = LFT3.Atom_Parsing_List(compl2_base, compl2_linker, compl1_base)
    id_dist_compl2 = LFT2.retrieve_atom_index(compl2_linker, APL1[2]) + compl1_base.mol_length
    id_dist_compl1 = LFT2.retrieve_atom_index(compl1_base, APL1[3])

    # Now we iterate over the array. When we calculate for a small distance, we will remember the index of its array and keep that as the best fit so far
    dist = 100      # start with an unreasonable distance. This value will be tested against the distances we find
    for i in range(len(compl1_base_confs)):
        for j in range(len(compl2_base_confs)):
            tmp_array = array_of_possible_dinucleotide_conformations[i][j]
            v_compl1 = tmp_array[id_dist_compl1]
            v_compl2 = tmp_array[id_dist_compl2]

            distance_between_nucleotides = LFT1.get_length_of_vector(v_compl1, v_compl2)

            if 1.25 <= distance_between_nucleotides < dist :
                # Override the variable dist and save the best dinucleotide conformation as an array. Also save off the index to parse the correct json file.
                dist = distance_between_nucleotides
                arr_best_fit = array_of_possible_dinucleotide_conformations[i][j]

                compl1_id = i
                compl2_id = j

    # If the distance is already of the desired size, then return the dinucleotide as is
    if LFT1.assert_length_of_vector(dist):
        complementary_strand = arr_best_fit
        return complementary_strand, index_lead

    # Now that we have the best fit, just to be sure we instantiate the correct object for the correct json filem for atom parsing.
    # This is generally not necessary, if the only differences between the json files, of the conformations of the same nucleosides, is the atoms coordinates.
    compl1_base = Nucleoside(compl1_base_confs[i])
    compl2_base = Nucleoside(compl2_base_confs[j])
    # get the size of the array as an index to correctly parse the required atoms of compl2_base
    index_fit = compl1_base.mol_length

    # Parse the two nucleotides from the array
    nuc1 = arr_best_fit[:index_fit]
    nuc2 = arr_best_fit[index_fit:]

    id_dist_compl2 = LFT2.retrieve_atom_index(compl2_linker, APL1[2])
    id_dist_compl1 = LFT2.retrieve_atom_index(compl1_base, APL1[3])
    ## GET THE VECTORS AND THE INDEX OF THE ATOMS REQUIRED FOR THE ROTATIONS
    # NUCLEOTIDE 1
    nuc1_base = compl1_base.get_base_denominator()
    nuc1_origin, _ = LFT3.retrieve_atoms_for_plane_rotation_of_complement(nuc1_base, nuc1_base)
    id_nuc1_origin = LFT2.retrieve_atom_index(compl1_base, nuc1_origin[0])

    nuc1_for_direction, _, __= LFT3.retrieve_atoms_for_positioning_of_complement1(nuc1_base, nuc1_base)
    id_nuc1_for_direction = LFT2.retrieve_atom_index(compl1_base, nuc1_for_direction[2])

    v_direction1 = LFT1.return_normalized(nuc1[id_nuc1_for_direction] - nuc1[id_nuc1_origin])
    # NUCLEOTIDE 2
    nuc2_base = compl2_base.get_base_denominator()
    nuc2_origin, _ = LFT3.retrieve_atoms_for_plane_rotation_of_complement(nuc2_base, nuc2_base)
    id_nuc2_origin = LFT2.retrieve_atom_index(compl2_base, nuc2_origin[0]) + compl2_linker.mol_length

    nuc2_for_direction, _ , __= LFT3.retrieve_atoms_for_positioning_of_complement1(nuc2_base, nuc2_base)
    id_nuc2_for_direction = LFT2.retrieve_atom_index(compl2_base, nuc2_for_direction[2]) + compl2_linker.mol_length

    v_direction2 = LFT1.return_normalized(nuc2[id_nuc2_for_direction] - nuc2[id_nuc2_origin])

    ## ASSERT THE ROTATIONS
    # use a small angle of 2 degrees
    angle_of_rot_TEST = 5 * (np.pi/180)

    # Rotate both directions the nuc1 and assert the distance between nuc1 and nuc2, then remember the direction of the smallest turn
    quat_nuc1 = LFT1.get_custom_quaternion(v_direction1, angle_of_rot_TEST)
    testnuc1 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc1, nuc1, nuc1[id_nuc1_origin])
    distance_between_nucs = LFT1.get_length_of_vector(testnuc1[id_dist_compl1], nuc2[id_dist_compl2])
    if not LFT1.assert_length_of_vector(distance_between_nucs) :
        x1 = distance_between_nucs

        v_direction1_test = v_direction1 * -1.0
        quat_nuc1 = LFT1.get_custom_quaternion(v_direction1_test, angle_of_rot_TEST)
        testnuc1 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc1, nuc1, nuc1[id_nuc1_origin])
        x2 = LFT1.get_length_of_vector(testnuc1[id_dist_compl1], nuc2[id_dist_compl2])
        if x1 > x2 :
            v_direction1 *= -1.0


    quat_nuc2 = LFT1.get_custom_quaternion(v_direction2, angle_of_rot_TEST)
    testnuc2 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc2, nuc2, nuc2[id_nuc2_origin])
    distance_between_nucs = LFT1.get_length_of_vector(nuc1[id_dist_compl1], testnuc2[id_dist_compl2])

    if not LFT1.assert_length_of_vector(distance_between_nucs) :
        x1 = distance_between_nucs

        v_direction2_test = v_direction2 * -1.0
        quat_nuc2 = LFT1.get_custom_quaternion(v_direction2_test, angle_of_rot_TEST)
        testnuc2 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc2, nuc2, nuc2[id_nuc2_origin])
        x2 = LFT1.get_length_of_vector(nuc1[id_dist_compl1], testnuc2[id_dist_compl2])
        if x1 > x2 :
            v_direction2 *= -1.0

    # Generate angles of rotation. Make note that we won't go past a 15 degree angle, since that general will not be necessary to turn that much.
    array_of_rot_angles = np.linspace(2.5, 30, 16) * (np.pi/180)

    stored_distances = np.zeros(len(array_of_rot_angles))
    for i in range(len(array_of_rot_angles)) :
        quat_nuc1 = LFT1.get_custom_quaternion(v_direction1, array_of_rot_angles[i])
        quat_nuc2 = LFT1.get_custom_quaternion(v_direction2, array_of_rot_angles[i])

        testnuc1 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc1, nuc1, nuc1[id_nuc1_origin])
        testnuc2 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc2, nuc2, nuc2[id_nuc2_origin])

        distance_between_nucs = LFT1.get_length_of_vector(testnuc1[id_dist_compl1], testnuc2[id_dist_compl2])
        # If the distance is suitable according to the boundaries, return the stack
        if LFT1.assert_length_of_vector(distance_between_nucs) :
            complementary_strand = np.vstack((testnuc1, testnuc2))
            return complementary_strand, index_lead
        else:
            stored_distances[i] = distance_between_nucs

    # If this part of the code is reached, that means none of the evaluated distances are suitable. Let's search the shortest one in the list
    min_val = stored_distances.min()
    index_min = np.where(stored_distances == min_val)
    angle_of_rot = array_of_rot_angles[index_min[0]]
    quat_nuc1 = LFT1.get_custom_quaternion(v_direction1, float(angle_of_rot))
    quat_nuc2 = LFT1.get_custom_quaternion(v_direction2, float(angle_of_rot))

    nuc1 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc1, nuc1, nuc1[id_nuc1_origin])
    nuc2 = LFT1.move_to_origin_ROTATE_move_back_to_loc(quat_nuc2, nuc2, nuc2[id_nuc2_origin])
    complementary_strand = np.vstack((nuc1, nuc2))
    return complementary_strand, index_lead


def assert_possible_base_conformations_and_fit(leading_nuc, leading_array : np.ndarray, conformations : list, compl_linker, complementary_strand : np.ndarray,
                                                                                prev_compl_nuc, prev_compl_linker, index_lead : int, index_compl : int) -> Union[np.ndarray, str]:
    """ Conformations is a list of the complementary_codex[NA], that contains one, two or three different conformations of the same NA.

        Conformations contains the names of particular json files, which will be converted to json objects.
        leading_nuc is already a json object.
        compl_linker is already a json object
        prev_compl_nuc is already json objects

        after the initial fit is done, it will calculate the best orientation for the nucleoside to be in, with respect to leading and complementary strand.
        Returns the conformation that has the best fit into the leading strand and the progressing complementary strand """

    # If there is only one conformation, there is no need to assert the differences in distance or dihedral
    if len(conformations) == 1:
        conf_n = Nucleoside(conformations[0])
        conf_n_arr = position_complementary_base(leading_nuc, conf_n, leading_array, index_lead)
        conf_with_linker = position_phosphate_linker(conf_n, conf_n_arr, compl_linker)
#        return conf_with_linker
        # Fit the nucleoside array better
        conf_n_final = tilt_array_to_get_a_better_fit(conf_n, compl_linker, prev_compl_nuc, prev_compl_linker, conf_with_linker, complementary_strand, index_compl)
        return conf_n_final

    ## Since there are multiple conformations available, we will need to sort out which one will fit the best.
    ## First we parse the correct indexes for the vectors we want to evaluate (check the distance)
    ## We first check if there are nucleotides that fit well without rotation of the base-plane. If that's not the case, we will check all conformations based on best fit after rotation of the base-plane.
    ## We then assert the different conformations after they've been fitted to the best it can

    # nuc_data is used just to be able to parse the correct indexes in the code below
    nuc_data = Nucleoside(conformations[0])

    # The first atom in the backbone of the nucleoside the current nucleoside attaches to
    atomOfInterest2 = LFT3.backbone_codex[json.loads(nuc_data.jason["identity"])[1]][0]   # first value of the backbone_codex to parse. else make list(_VAL)[0]
    bb_id = LFT2.retrieve_atom_index(prev_compl_nuc, atomOfInterest2) + index_compl + prev_compl_linker.mol_length
    bb_v = complementary_strand[bb_id]

    # The last atom in the backbone of the linker of the current nucleoside. 
    atomOfInterest1 = LFT3.backbone_codex[json.loads(compl_linker.jason["identity"])[0]][-1] # since it is the last of the backbone_codex that should be parsed. else make list(_VAL)[-1]
    link_id = LFT2.retrieve_atom_index(compl_linker, atomOfInterest1)

    ## First we assert the conformations without tilting by the base-plane

    # Instantiate a matrix to store possible conformations in an array
    possibilities_of_conformations = np.zeros(shape=len(conformations), dtype=object)
    # Instantiate a matrix to store the distances in the backbone, from linker to previous nucleotide, and whether or not they fit within the boundaries of asserting set distances
    stored_bb_distances = np.zeros(shape=len(conformations), dtype=object)
    stored_bb_distance_bools = np.zeros(shape=len(conformations), dtype=bool)
    stored_bb_dihedral_bools = np.zeros(shape=len(conformations), dtype=bool)


    for file_n in range(len(possibilities_of_conformations)):
        conf_n = Nucleoside(conformations[file_n])
        conf_n_arr = position_complementary_base(leading_nuc, conf_n, leading_array, index_lead)
        # add the linker to the nucleoside
        possibilities_of_conformations[file_n] = position_phosphate_linker(conf_n, conf_n_arr, compl_linker)

        link_v = possibilities_of_conformations[file_n][link_id]
        # Check the distance between the two nucleotides, so P -> O3' (with native nucs as example) one of which will have the shortest distance. Pick that one
        stored_bb_distances[file_n] = LFT1.get_length_of_vector(link_v, bb_v)
        stored_bb_distance_bools[file_n] = LFT1.assert_length_of_vector(stored_bb_distances[file_n])
        stored_bb_dihedral_bools[file_n] = LFFB.assert_the_dihedral_of_interest(conf_n, possibilities_of_conformations[file_n], compl_linker, prev_compl_nuc, complementary_strand, index_compl, prev_compl_linker)

    # check dihedral suitability
    if np.any(stored_bb_dihedral_bools == True) :
        index_of_best_conformation = LFFB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_dihedral_bools)
        if not stored_bb_distance_bools[index_of_best_conformation] == True:
            return tilt_array_to_get_a_better_fit(Nucleoside(conformations[index_of_best_conformation]), compl_linker, prev_compl_nuc, prev_compl_linker, possibilities_of_conformations[index_of_best_conformation], complementary_strand, index_compl)

        return possibilities_of_conformations[index_of_best_conformation]
    # if the check dihedral suitability fails, check for distance suitability
    if np.any(stored_bb_distance_bools == True):
        index_of_best_conformation = LFFB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_distance_bools)
        return possibilities_of_conformations[index_of_best_conformation]


    ## Alas, none of the conformations were suitable enough, it seems we will have to fit them over and over.
    # Instantiate a matrix to store possible conformations in an array
    possibilities_of_conformations = np.zeros(shape=len(conformations), dtype=object)
    # Instantiate a matrix to store the distances in the backbone, from linker to previous nucleotide, and whether or not they fit within the boundaries of asserting set distances
    stored_bb_distances = np.zeros(shape=len(conformations), dtype=object)
    stored_bb_distance_bools = np.zeros(shape=len(conformations), dtype=bool)
    stored_bb_dihedral_bools = np.zeros(shape=len(conformations), dtype=bool)

    for file_n in range(len(possibilities_of_conformations)):
        conf_n = Nucleoside(conformations[file_n])
        conf_n_arr = position_complementary_base(leading_nuc, conf_n, leading_array, index_lead)
        # add the linker to the nucleoside
        conf_n_arr_link = position_phosphate_linker(conf_n, conf_n_arr, compl_linker)
        # tilt the conformation to the best possible fit
        possibilities_of_conformations[file_n] = tilt_array_to_get_a_better_fit(conf_n, compl_linker, prev_compl_nuc, prev_compl_linker, conf_n_arr_link, complementary_strand, index_compl)

        link_v = possibilities_of_conformations[file_n][link_id]
        # Check the distance between the two nucleotides, so P -> O3' (with native nucs as example) one of which will have the shortest distance. Pick that one
        stored_bb_distances[file_n] = LFT1.get_length_of_vector(link_v, bb_v)
        stored_bb_distance_bools[file_n] = LFT1.assert_length_of_vector(stored_bb_distances[file_n])
        stored_bb_dihedral_bools[file_n] = LFFB.assert_the_dihedral_of_interest(conf_n, possibilities_of_conformations[file_n], compl_linker, prev_compl_nuc, complementary_strand, index_compl, prev_compl_linker)

    if np.any(stored_bb_dihedral_bools == True) :
        index_of_best_conformation = LFFB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_dihedral_bools)
        return possibilities_of_conformations[index_of_best_conformation]
    # This gives the one with the least distance from the desired bb_distance
    index_of_best_conformation = LFFB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_distance_bools)
    return possibilities_of_conformations[index_of_best_conformation]


def generate_complementary_sequence(sequence_list : list, complement : Union[list, str]) -> list:
    """ sequence list is the given input.
        complement will specify what the complementary strand will look like. choices between homo - DNA - RNA """
    complementary_dictDNA = { "A" : "T", "T" : "A", "G" : "C", "C" : "G", "U" : "A" }
    complementary_dictRNA = { "A" : "U", "T" : "A", "G" : "C", "C" : "G", "U" : "A" }

    bases = LFT2.retrieve_bases_list(sequence_list)
    if isinstance(complement, list):
        import fundaments

        complementary_sequence = list(map(lambda x: x.strip(","), complement))

        # See of the lengths match. If the lengths do not match, give an assertion error and print the following string. This exits the software too.
        assert len(sequence_list) == len(complementary_sequence), "The length of the complementary strand does not match the length of the leading strand!"

        # Check if one of the nucleosides in the prompted list is wrong, i.e. not existing or wrongly prompted (misspelled)

        # The keys, meaning the nucleosides, from the complementary dictionary wil be parsed as a list
        keys_of_dict = LFT3.codex_acidum_nucleicum.keys()
        if not fundaments.check_if_nucleotides_are_valid(keys_of_dict):
            sys.exit(0)

        # Reverse the prompted list
        return complementary_sequence[::-1]

    # if the sequence needs to be complementary per nucleoside
    if complement.lower() == "homo":
        complementary_sequence = []
        list_of_chemistries = LFT2.retrieve_chemistry_list(sequence_list)

        # The keys, meaning the nucleosides, from the complementary dictionary wil be parsed as a list
        keys_of_dict = LFT3.codex_acidum_nucleicum.keys()

        # Get only the nucleosides of the same chemistry as the chemistry it will complement
        for chem_i in range(len(list_of_chemistries)):
            # Get the possible nucleosides we can choose from with that have the same chemistry as in the leadstrand's nucleoside
            parsed_chemistry, ln_str = LFT2.retrieve_chemistry(list_of_chemistries[chem_i])
            homo_chemistry_list = LFT2.retrieve_homo_nucleosides(keys_of_dict, parsed_chemistry, ln_str)

            # Assess which base fits the complementary nucleoside
            complementary_nucleoside = LFT2.assess_possible_complementary_base(parsed_chemistry, bases[chem_i], homo_chemistry_list, complementary_dictRNA, complementary_dictDNA)
            complementary_sequence.append(complementary_nucleoside)

        return complementary_sequence

        # We will leave this part of the old code here just in case, but I don't think we'll need it anymore
            #complementary_sequence = []
            #for chem_i in range(len(chemistry)):
            #    if chemistry[chem_i] == "r":
            #        comp_base = complementary_dictRNA[bases[chem_i]]
            #        complementary_sequence.append(chemistry[chem_i] + comp_base)
            #    else:
            #        comp_base = complementary_dictDNA[bases[chem_i]]
            #        complementary_sequence.append(chemistry[chem_i] + comp_base)

            #return complementary_sequence


    if complement.upper() == "DNA":
        chemistry = "d"

        # Switch the bases the get their complementary base
        comp_bases = LFT2.get_complementary_bases(bases, complementary_dictDNA)
        complementary_sequence = LFT2.concatenate_chem_and_bases(chemistry, comp_bases)
        return complementary_sequence

    if complement.upper() == "RNA":
        chemistry = "r"

        # Switch the bases the get their complementary base
        comp_bases = LFT2.get_complementary_bases(bases, complementary_dictRNA)
        complementary_sequence = LFT2.concatenate_chem_and_bases(chemistry, comp_bases)
        return complementary_sequence

    if True:
        # At this point, any input the user has prompted should have gone through a return statement. So if we reach this point, just stop the program. 
        raise ValueError("The variable you have prompted for the '--complement' flag is not correct. Please review your input file.\n")
        sys.exit(1)


def cap_nucleic_acid_strands(leading_array : np.ndarray, leading_sequence : list, complementary_array : np.ndarray, complementary_sequence : list) -> Union[np.ndarray, list]:
    """ Cap the nucleic acid strands with a hydrogen, to finish the build of the duplex
        We create two functions in LFFB that parse both the correct coordinates of the capping atoms and their names """

    # Retrieve the backbone dictionary from LFT3, since we'll be needing it
    backbone = LFT3.backbone_codex

    # Retrieve the dictionary for the nucleoside filenames from LFT3 as well
    filenames = LFT3.codex_acidum_nucleicum

    # Defines and returns the cartesian position of the hydrogens 
    atom_array = LFFB.capping_retrieve_atomarrays(leading_array, leading_sequence, complementary_array, complementary_sequence, backbone, filenames)

    # Returns the name of the hydrogens that have been defined with the previous function
    atom_names = LFFB.capping_retrieve_atomnames(leading_sequence, complementary_sequence, backbone, filenames)

    return atom_array, atom_names


def create_PDB_from_array_final(leading_array : np.ndarray, list_of_leading_sequence : list, complementary_array : np.ndarray, list_of_complementary_sequence : list):
    """ Write out the data for the pdb filename
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html"""

    # Reverse the sequence of this to parse the correct atoms along the array
    list_of_leading_sequence = list_of_leading_sequence[::-1]

    # Capping of the duplexes
    # First two elements belong to the leading strand, the last two elements belong to the complementary strand. A total of four elements per variable
    capping_array, capping_names = cap_nucleic_acid_strands(leading_array, list_of_leading_sequence, complementary_array, list_of_complementary_sequence)

    # Adjust leading strand's array for the capping
    leading_array = np.vstack((capping_array[0], leading_array, capping_array[1]))

    # LEADING STRAND
    df_leading = pd.DataFrame()

    df_leading["RecName"] = ["ATOM" for x in range(leading_array.shape[0])]
    df_leading["AtomNum"] = np.arange(start=1, stop=leading_array.shape[0] + 1)
    df_leading["AtomName"] = capping_names[0] + LFT2.LEAD_pdb_AtomNames_or_ElementSymbol(list_of_leading_sequence, "Atoms") + capping_names[1]
    df_leading["AltLoc"] = " "
    df_leading["ResName"] = LFT2.LEAD_pdb_Residuename(list_of_leading_sequence)
    df_leading["Chain"] = "J"
    df_leading["Sequence"] = LFT2.LEAD_pdb_Sequence(list_of_leading_sequence)
    df_leading["X_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,0]))
    df_leading["Y_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,1]))
    df_leading["Z_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,2]))
    df_leading["Occupancy"] = "1.00"
    df_leading["Temp"] = "0.00"
    df_leading["SegmentID"] = str("   ")
    df_leading["ElementSymbol"] = ["H"] + LFT2.LEAD_pdb_AtomNames_or_ElementSymbol(list_of_leading_sequence, "Symbol") + ["H"]

    # TER line between the single strands
    TER_line = pd.DataFrame({"RecName" : "TER", "AtomNum" : 69}, index=[69])
    ln_lead = len(df_leading["AtomNum"])

    # Adjust complementary strand's array for the capping
    complementary_array = np.vstack((capping_array[2], complementary_array, capping_array[3]))

    # COMPLEMENTARY STRAND
    df_complementary = pd.DataFrame()

    df_complementary["RecName"] = ["ATOM" for x in range(complementary_array.shape[0])]
    df_complementary["AtomNum"] = np.arange(start=ln_lead + 1, stop=(ln_lead + complementary_array.shape[0] + 1))
    df_complementary["AtomName"] = capping_names[2] + LFT2.COMPLEMENTARY_pdb_AtomNames_or_ElementSymbol(list_of_complementary_sequence, "Atoms") + capping_names[3]
    df_complementary["AltLoc"] = " "
    df_complementary["ResName"] = LFT2.COMPLEMENTARY_pdb_Residuename(list_of_complementary_sequence)
    df_complementary["Chain"] = "R"
    df_complementary["Sequence"] = LFT2.COMPLEMENTARY_pdb_Sequence(list_of_complementary_sequence, len(list_of_complementary_sequence))
    df_complementary["X_coord"] = list(map(lambda x: "{:.3f}".format(x), complementary_array[:,0]))
    df_complementary["Y_coord"] = list(map(lambda x: "{:.3f}".format(x), complementary_array[:,1]))
    df_complementary["Z_coord"] = list(map(lambda x: "{:.3f}".format(x), complementary_array[:,2]))
    df_complementary["Occupancy"] = "1.00"
    df_complementary["Temp"] = "0.00"
    df_complementary["SegmentID"] = str("   ")
    df_complementary["ElementSymbol"] = ["H"] + LFT2.COMPLEMENTARY_pdb_AtomNames_or_ElementSymbol(list_of_complementary_sequence, "Symbol") + ["H"]

    # Concatenate the two dataframes
    duplex_df = pd.concat([df_leading, TER_line, df_complementary], ignore_index=True)

    # Write out the pdb file
    filename = "testing_duplex.pdb"
    with open(filename ,"w") as pdb:
        for index, row in duplex_df.iterrows():
            if row[0] == "TER":
                TERline = row[0]
                pdb.write("%-6s\n" % TERline)
            if not row[0] == "TER":
                split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
                pdb.write("%-6s%5s%5s%s%3s%2s%5d  %8s%8s%9s%6s%7s%4s     %2s\n" % tuple(split_line))
        pdb.write("END")
        pdb.close()

    print("\n\nWriting nucleic acid duplex to " + filename + ". \n\n")



#def create_PDB_from_array(leading_array : np.ndarray, list_of_leading_sequence : list) -> None:
#    """ Write out the data for the pdb file """
#    print("Writing to pdb ...")
#    list_of_leading_sequence = list_of_leading_sequence[::-1]
#
#    df = pd.DataFrame()
#
#    df["RecName"] = ["ATOM" for x in range(leading_array.shape[0])]
#    df["AtomNum"] = np.arange(start=1, stop=leading_array.shape[0] + 1)
#    df["AtomName"] = LFT2.pdb_AtomNames_or_ElementSymbol(list_of_leading_sequence, "Atoms")
#    df["AltLoc"] = " "
#    df["ResName"] = LFT2.pdb_Residuename(list_of_leading_sequence)
#    df["Chain"] = "A"
#    df["Sequence"] = LFT2.pdb_Sequence(list_of_leading_sequence)
#    df["X_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,0]))
#    df["Y_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,1]))
#    df["Z_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,2]))
#    df["Occupancy"] = "1.00"
#    df["Temp"] = "0.00"
#    df["SegmentID"] = str("   ")
#    df["ElementSymbol"] = LFT2.pdb_AtomNames_or_ElementSymbol(list_of_leading_sequence, "Symbol")
#
#
#    # Write out the pdb file
#    filename = "testing_daedalus.pdb"
#    with open(filename ,"w") as pdb:
#        for index, row in df_leading.iterrows():
#            split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
#            pdb.write("%-6s%5s%5s%s%2s%3s%5s  %8s%8s%8s%6s%6s%4s      %2s\n" % tuple(split_line))
#        pdb.write("END")
#        pdb.close()
