import numpy as np
from typing import Union

import initMolecule

import builder.utils_builder as UB
import builder.mathematics as MATH
import builder.parse_or_write as PARSE
import builder.builder_constants as CONSTANTS



""" This file is used to refactor some of the code from utils_builder.py, to make a more legible. """

def retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_bools):

    # If this case happens, then check which one is closest to the defaulted size between the linker and the previous nucleotide
    if np.all(stored_bb_bools == True):
        diff_distances = np.zeros(shape=len(stored_bb_bools), dtype=float)
        for i in range(len(diff_distances)):
            diff_distances[i] = abs(stored_bb_distances[i] - 1.6)

        smallest_diff = diff_distances.min()
        index_smallest_dif = np.where(diff_distances == smallest_diff)[0]
        return int(index_smallest_dif)

    # Check if some are True. If there is only one True, just take that conformation. If there are multiple Trues, test them against eachother
    if np.any(stored_bb_bools == True):
        if np.count_nonzero(stored_bb_bools == True) == 1:
            only_true_conf = np.where(True == stored_bb_bools)[0]
            return int(only_true_conf)

        multiple_true_conf = np.where(True == stored_bb_distances)[0]
        diff_distances = np.zeros(shape=len(multiple_true_conf), dtype=float)
        for i in multiple_true_conf:
            diff_distances[i] = abs(stored_bb_distances[i] - 1.6)

        smallest_diff = diff_distances.min()
        index_smallest_dif = np.where(diff_distances == smallest_diff)[0]
        return int(index_smallest_dif)

    # if the program reaches this part, no distances were suitable. Just test them all against the 1.6 and get the smallest value of it, return that
    diff_distances = np.zeros(shape=len(stored_bb_bools), dtype=float)
    for i in range(len(diff_distances)):
        diff_distances[i] = abs(stored_bb_distances[i] - 1.6)

    smallest_diff = diff_distances.min()
    index_smallest_dif = np.where(diff_distances == smallest_diff)[0]
    return int(index_smallest_dif)


def return_cross_vector_for_plane_rotation(strand_array : np.ndarray, atom_list_id : list) -> np.ndarray :
    """ Take the array of either the leading or complementary strand and the indexes of the atoms we want to retrieve from it.
        returns the cross product of the of the plane in which the base of the nucleoside lies in.
        In other words, it give you a vector that is perpendicular to the nucleoside's base. """

    v0 = strand_array[atom_list_id[0]]      # Atom that connects base to the sugar moiety
    v1 = strand_array[atom_list_id[1]]      # Atom that, with the next vector v2, defines the plane in which the base lies in.
    v2 = strand_array[atom_list_id[2]]      # Atom that, with the previous vector v1, defines the plane in which the base lies in.

    p0 = MATH.return_normalized(v1 - v0)    # Vector connect -> v1
    p1 = MATH.return_normalized(v2 - v0)    # Vector connect -> v2

    # Return the cross vector of the plane in which the nucleoside's lies in.
    return MATH.return_normalized(np.cross(p1, p0))


def return_position_for_complementary_base(angle : float, dihedral : float, strand_array : np.ndarray, atom_list_id : list, vector_size : float) -> np.ndarray : 
    """ Returns the vector that denotes the position of the where we want to complementary base to rotate to"""

    # Get the vectors of interest from the array
    v0 = strand_array[atom_list_id[0]]
    v1 = strand_array[atom_list_id[1]]
    v2 = strand_array[atom_list_id[2]]

    # Generate a vector of interest, here in the context of positioning the complementary nucleoside's base
    v_position_for_complement = UB.generate_vector_of_interest(angle, dihedral, [v2, v1, v0])
    return MATH.return_normalized(v_position_for_complement) * vector_size


def position_next_nucleoside(next_nucleoside, prev_nucleoside, prev_linker, leading_strand : np.ndarray) -> np.ndarray:
    """ This function is used after position_phosphate_linker().
        It serves the purpose of adding the next nucleotide onto the leading_strand.

        next_nucleoside : json object
        leading_strand : the nucleic acid strand to which we append next_nucleoside to. """
    # The first thing to do is to find the location of the subsequent atom, here O3', then rotate the next_nucleoside by zeta and epsilon
    # Afterwards, we turn over the epsilon dihedral and by rotate the normal of the plane have to the plane we want; this positions everything!

    # Atom Parsing List (ATP) = Parse which linker and which nucleotide the previous one is
    APL = PARSE.Atom_Parsing_List(prev_nucleoside, prev_linker, next_nucleoside)
    # Dihedral Parsing List (DPL) = Parse which dihedrals are required to rotate on and over
    # Angle Parsing List (AngPL)
    AngPL = DPL =  PARSE.retrieve_list_of_dihedrals_and_angles_to_build_with(next_nucleoside)

    #### POSITION THE NEXT NUCLEOSIDE PROPERLY.
    # dihedral C5' - O5' - P - O3'
    # get the size of the linker, since the array goes : " linker.array + nucleotide.array", where linker.array precedes the indices of the nucleotide.array
    shape_linker = prev_linker.mol_length
    id_v0 = PARSE.retrieve_atom_index(prev_nucleoside, APL[0]) + shape_linker
    id_v1 = PARSE.retrieve_atom_index(prev_nucleoside, APL[1]) + shape_linker
    id_v2 = PARSE.retrieve_atom_index(prev_linker,APL[2])
    v0 = leading_strand[id_v0]
    v1 = leading_strand[id_v1]
    v2 = leading_strand[id_v2]

    # Get angle and dihedral:
    alpha_angle = next_nucleoside.get_angle(AngPL[0])
    alpha_dihr = next_nucleoside.get_dihedral(DPL[0])

    ## We position the nextnuc by the position of O3'
    single_vector1 = UB.generate_vector_of_interest(alpha_angle, alpha_dihr, [v2, v1, v0])
    linker_nuc_distance = 1.6
    # Add single_vector1 to v2, so that we define the location of ATP[3]
    v3 = MATH.move_vector_to_loc(MATH.return_normalized(single_vector1) * linker_nuc_distance, v2)

    ## Bring the nucleoside array to the new position
    id_v3 = PARSE.retrieve_atom_index(next_nucleoside, APL[3])
    v3_original = next_nucleoside.array[id_v3]
    # Get distance from next_nucleoside O3' (v3_old) to the position defined as O3' (v3), so the distance that we need to move the array with
    distance_v3_v3_original = v3 - v3_original
    # Get the distance from the P -> O3' now and add it to nextnuc_origin which will move the entire molecule(nextnuc) to the position of O3'
    next_nucleoside_loc = MATH.move_vector_to_loc(next_nucleoside.array, distance_v3_v3_original)

    for i in range(1, len(DPL)):
        # Only at the final dihedral do we require a plane rotation.
        # The reason we forloop is because perhaps there might be cases where just a single rotation and a plane rotation are sufficient, like with morpholino's presumably.
        if i + 1 == len(DPL):
            # Do the final vector rotation + plane rotation
            id_vA = PARSE.retrieve_atom_index(prev_linker, APL[i])
            id_vB = PARSE.retrieve_atom_index(next_nucleoside, APL[i + 1])
            id_vC = PARSE.retrieve_atom_index(next_nucleoside, APL[i + 2])
            vA = leading_strand[id_vA]
            vB = next_nucleoside_loc[id_vB]
            vC = next_nucleoside_loc[id_vC]
            #Get the required angles
            angle_N = next_nucleoside.get_angle(AngPL[i])
            dihedral_N = next_nucleoside.get_dihedral(DPL[i])

            single_vector_N = UB.generate_vector_of_interest(angle_N, dihedral_N, [vC, vB, vA])

            # Retrieve the next atom vector in the sequence, required for rotation
            id_vD = PARSE.retrieve_atom_index(next_nucleoside, APL[i + 3])
            vD = next_nucleoside_loc[id_vD]

            ## Generate to normal vectors of the planes of interest. The order in which you perform the cross product is not important BUT!!!
            # It IS important that the same vectors are operated on in the same order for both normal vectors!!!
            # normal to rotate from
            n0 = MATH.get_normal_vector_of_plane(vC - vB, vD - vC)
            # normal to rotate to
            n1 = MATH.get_normal_vector_of_plane(vC - vB, single_vector_N)
            # The order of the quaternion does matter, as it starts with " vector to rotate to" and secondly with "vector that we want to rotate from "
            quaternion_plane = MATH.get_quaternion(n1, n0)

            #Now bring the array to the origin
            distance_to_origin_N = next_nucleoside_loc[id_vC]
            next_nucleoside_loc = MATH.move_to_origin_ROTATE_move_back_to_loc(quaternion_plane, next_nucleoside_loc, distance_to_origin_N)

            return next_nucleoside_loc

        # Just continue if it is not the last one in the loop
        id_vA = PARSE.retrieve_atom_index(prev_nucleoside, APL[i]) + shape_linker
        vA = leading_strand[id_vA]
        id_vB = PARSE.retrieve_atom_index(prev_linker, APL[i + 1])
        vB = leading_strand[id_vB]
        id_vC = PARSE.retrieve_atom_index(next_nucleoside, APL[i + 2])
        vC = next_nucleoside_loc[id_vC]

        #Get the required angles
        angle_N = next_nucleoside.get_angle(AngPL[i])
        dihedral_N = next_nucleoside.get_dihedral(DPL[i])

        single_vector_N = UB.generate_vector_of_interest(angle_N, dihedral_N, [vC, vB, vA])

        # Retrieve the appropriate quaternion for the rotation 
        id_vD = PARSE.retrieve_atom_index(next_nucleoside, APL[i + 3])
        vD = next_nucleoside_loc[id_vD]
        pC_D = MATH.return_normalized(vD - vC)
        quaternion_N = MATH.get_quaternion(single_vector_N, pC_D)

        #Now rotate your molecule onto single_vector_N
        distance_to_origin_N = next_nucleoside_loc[id_vC]

        next_nucleoside_loc = MATH.move_to_origin_ROTATE_move_back_to_loc(quaternion_N, next_nucleoside_loc, distance_to_origin_N)


def assert_the_dihedral_of_interest(compl_nuc, compl_nuc_arr : np.ndarray, compl_linker, prev_nuc, complementary_strand : np.ndarray, index_compl, prev_linker) -> bool:
    # Atom Parsing List for the knowing which atoms to parse from the respective arrays ; for bond length and dihedral evaluation
    APL = PARSE.Atom_Parsing_List(compl_nuc, compl_linker, prev_nuc)
    dihedral = compl_nuc.get_dihedral("alpha")
    angle = compl_nuc.get_angle("alpha")

    # Parse the last atom needed from the complementary strand. NB : index_compl has been set to the correct value in the function ' assert_possible ... _and_fit() '.
    id_compl_strand = PARSE.retrieve_atom_index(prev_nuc, APL[3]) + index_compl + prev_linker.mol_length
    v_compl_strand = complementary_strand[id_compl_strand]

    # Parse the other atoms needed from the current nucleotide that is being fitted
    id_v0 = PARSE.retrieve_atom_index(compl_nuc, APL[0]) + compl_linker.mol_length
    id_v1 = PARSE.retrieve_atom_index(compl_nuc, APL[1]) + compl_linker.mol_length
    id_v2 = PARSE.retrieve_atom_index(compl_linker, APL[2])

    v0 = compl_nuc_arr[id_v0]
    v1 = compl_nuc_arr[id_v1]
    v2 = compl_nuc_arr[id_v2]

    # Calculate the dihedral. Let's say the dihedral should not deviate more than 25 degrees?
    calculated_dihr = MATH.dihedral_single(v0, v1, v2, v_compl_strand)

    # Create booleans to assert the calculated values and see if they are within ranges of concordance
    bool_dihedral = MATH.assert_dihedral_of_nucleotide(calculated_dihr, dihedral)

    return bool_dihedral


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
    rotPlane_leadingBase_atoms, rotPlane_complBase_atoms = CONSTANTS.retrieve_atoms_for_plane_rotation_of_complement(leadingBase, complBase)
    ## Get the proper vectors for the plane rotation
    # Leading base vector
    rotPlane_leadingBase_atoms_id = PARSE.retrieve_atom_index_MULTIPLE(leading_base, rotPlane_leadingBase_atoms, index_lead)
    vCross_leadingBase = return_cross_vector_for_plane_rotation(leading_array, rotPlane_leadingBase_atoms_id)

    # Complementary base vector
    rotPlane_complBase_atoms_id = PARSE.retrieve_atom_index_MULTIPLE(complementary_base, rotPlane_complBase_atoms)
    vCross_complBase = return_cross_vector_for_plane_rotation(complementary_base.array, rotPlane_complBase_atoms_id)

    ### ROTATE THE PLANE
    ## Rotate the complementary base onto the plane of the leading base
    # Get the cross product of the atoms of leading_base that make it if both bases are in the same plane, the cross product vectors are exactly opposite
    plane1_quaternion = MATH.get_quaternion(vCross_leadingBase * -1.0, vCross_complBase)

    ## Apply the rotation
    # vC0 is atom that binds the base with the sugar moiety, making it the central atom over which we rotate when we bring the molecule to the origin
    vC0 = complementary_base.array[rotPlane_complBase_atoms_id[0]]
    tmp_compl_arr = MATH.move_vector_to_origin(complementary_base.array, vC0)
    compl_nucleoside_array = MATH.rotate_with_quaternion(plane1_quaternion, tmp_compl_arr)

    # Get the required parameters for the rotations
    Q_angle, Q_dihedral, R_angle, R_dihedral = CONSTANTS.retrieve_angles_and_dihedrals_for_initial_base_positioning(leadingBase)

    ## Apply a translation to get the initial positioning
    Q_leadingBase_atoms, Q_complBase_atom, Q_distance = CONSTANTS.retrieve_atoms_for_positioning_of_complement1(leadingBase, complBase)
    Q_leadingBase_atoms_id = PARSE.retrieve_atom_index_MULTIPLE(leading_base, Q_leadingBase_atoms, index_lead)

    # Calculate the position at which we want to the R positioning to be. See paper @Figure X to show what R positioning is.
    Q_position = return_position_for_complementary_base(Q_angle, Q_dihedral, leading_array, Q_leadingBase_atoms_id, Q_distance)

    ## Now that we have the Q_position, we move the complementary base to the position by doing the following :
    # Get distance at which we want to position the base's hydrogen bond
    Q_vector = leading_array[Q_leadingBase_atoms_id[2]] + Q_position
    # retrieve the vector of the atom that we want to attach 'Q_vector' to.
    Q_complBase_atom_id = PARSE.retrieve_atom_index(complementary_base, Q_complBase_atom)

    # the distance from the location we want to have compl1_base2_atom move to. But we move the array from where it is at now, so substract that and then move it
    move_to_Q = Q_vector - compl_nucleoside_array[Q_complBase_atom_id]

    compl_nucleoside_array = MATH.move_vector_to_loc(compl_nucleoside_array, move_to_Q)

    ### After positioning the planes and moving the nucleoside to Q, we now rotate the nucleoside's base to R
    R_leadingBase_atoms, R_complBase_atom, R_distance = CONSTANTS.retrieve_atoms_for_position_of_complement2(leadingBase, complBase)
    R_leadingBase_atoms_id = PARSE.retrieve_atom_index_MULTIPLE(leading_base, R_leadingBase_atoms, index_lead)

    R_position = return_position_for_complementary_base(R_angle, R_dihedral, leading_array, R_leadingBase_atoms_id, R_distance)

    # Creates the coordinate at which we want to move our the Q_complement atom, from the complementary nucleoside's base to
    R_move_to = leading_array[R_leadingBase_atoms_id[2]] + R_position
    R_compl_id = PARSE.retrieve_atom_index(complementary_base, R_complBase_atom)
    R_vector = compl_nucleoside_array[Q_complBase_atom_id]

    # Create quaternion to rotate the complementary base a final time
    j1 = MATH.return_normalized(R_move_to - R_vector)
    j2 = MATH.return_normalized(compl_nucleoside_array[R_compl_id] - R_vector)

    # Check if the nucleoside needs reorientation. Apparently, if the vectors are nearing a dot product of -1, the rotations are pretty sketchy
    if MATH.assert_dot_product(np.dot(j1,j2)) :
        complementary_array = MATH.reorient_nucleoside_array(complementary_base.array)
        # Get first rotation
        vCross_complBase = return_cross_vector_for_plane_rotation(complementary_array, rotPlane_complBase_atoms_id)
        plane1_quaternion = MATH.get_quaternion(vCross_leadingBase * -1.0, vCross_complBase)
        # Apply rotation
        vC0 = complementary_array[rotPlane_complBase_atoms_id[0]]
        tmp_compl_arr = MATH.move_vector_to_origin(complementary_array, vC0)
        compl_nucleoside_array = MATH.rotate_with_quaternion(plane1_quaternion, tmp_compl_arr)
        # Denote location to move to
        move_to_Q = Q_vector - compl_nucleoside_array[Q_complBase_atom_id]
        compl_nucleoside_array = MATH.move_vector_to_loc(compl_nucleoside_array, move_to_Q)

        # Creates the coordinate at which we want to move our the Q_complement atom, from the complementary nucleoside's base to
        R_move_to = leading_array[R_leadingBase_atoms_id[2]] + R_position
        R_compl_id = PARSE.retrieve_atom_index(complementary_base, R_complBase_atom)
        R_vector = compl_nucleoside_array[Q_complBase_atom_id]

        # Create quaternion to rotate the complementary base a final time
        j1 = MATH.return_normalized(R_move_to - R_vector)
        j2 = MATH.return_normalized(compl_nucleoside_array[R_compl_id] - R_vector)

    # Final rotation and return the positioned array
    plane2_quaternion = MATH.get_quaternion(j1, j2)

    return MATH.move_to_origin_ROTATE_move_back_to_loc(plane2_quaternion, compl_nucleoside_array, R_vector)


def assess_complX_id(compl1_id : int, compl2_id : int) -> bool :
    """ If the defaulted values of the the variables remain as `-1`, that means the function was not passed a different value.
        This means that the variable complX_id was never assigned a proper index value for the array to begin with. This should not happen.

        Instead, we short circuit the Ducque program and ask the user to input different input parameters of the dihedrals, to see if the fits better."""

    try :
        if compl1_id == -1 or compl2_id == -1 :
            raise UnboundLocalError
    except UnboundLocalError :
            return False

    return True


def assert_rotation_of_bases_by_distance(array_nucs : list, v_directions: list, index_origin : Union[list, float], idxDistanceBetweenSubsequentNucs : list) -> np.ndarray:
    """ Assert the direction of the rotation axis by measuring the distances between the two nucleic acids; by the backbone"""
    ## ASSERT THE ROTATIONS
    nuc1 = array_nucs[0]
    nuc2 = array_nucs[1]

    # Initialise the values which we use for the ranges of rotation angles over which we rotate over
    _rot0 = np.deg2rad(5)
    _rotI = np.deg2rad(355)
    _rot00, _rot01, _rotI0, _rotI1 = np.deg2rad(2.5), np.deg2rad(20), np.deg2rad(357.5), np.deg2rad(340)

    # One set of nucleobase-plane rotations to check
    if len(v_directions) == 1:
        angles_of_rotation = np.array([_rot0, _rotI])
        v_direction = v_directions[0]
        #direction_of_rotation = np.array([[1], [-1]])
        _storedDistances = np.zeros(2)

        _returnAnglesOfRotation = np.array([np.linspace(_rot00, _rot01, 16),
                                           np.linspace(_rotI0, _rotI1, 16)])

        for i, angle in enumerate(angles_of_rotation) :
            # The direction of the direction_axis is the following ;
            #d1 = v_direction * set_of_angles[index]

            # Get quaternion for the rotation
            _quaternion = MATH.get_quaternion_custom_axis(v_direction, angle)

            # Rotate the nucleoside by the plane of the nucleobase
            _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternion, nuc2, nuc2[index_origin])

            # Calculate the rotation
            _storedDistances[i] = MATH.get_length_of_vector(nuc1[idxDistanceBetweenSubsequentNucs[1]], _rotatedNuc2[idxDistanceBetweenSubsequentNucs[0]])


        _storedDistances = MATH.smallest_difference(_storedDistances, 1.6)
        _indexMin = np.where(_storedDistances == _storedDistances.min())[0][0]
        # index_min == 0 is rotation 2.5 -> 30 degrees 
        # index_min == 1 is rotation 357.5 -> 330 degrees
        return _returnAnglesOfRotation[0] if _indexMin == 0 else _returnAnglesOfRotation[1]

        #direction = direction_of_rotation[index_min[0]][0]
        #return v_direction1 * direction

    # Two sets of nucleobase-plane rotations to check, when fitting the first two complementary nucleosides
    if len(v_directions) == 2:
        v_direction1, v_direction2 = v_directions[0], v_directions[1]
        _storedDistances = np.zeros(4)
        #direction_of_rotation = np.array([[1,1], [1,-1], [-1, 1], [-1, -1]])
        _anglesOfRotation = np.array([_rot0, _rot0,
                                      _rot0, _rotI,
                                      _rotI, _rot0,
                                      _rotI, _rotI])

        _returnAnglesOfRotation = np.array([np.linspace(_rot00, _rot01, 16), np.linspace(_rot00, _rot01, 16),
                                            np.linspace(_rot00, _rot01, 16), np.linspace(_rotI0, _rotI1, 16),
                                            np.linspace(_rotI0, _rotI1, 16), np.linspace(_rot00, _rot01, 16),
                                            np.linspace(_rotI0, _rotI1, 16), np.linspace(_rotI0, _rotI1, 16)])

        # Reshape the arrays to the desired (4,2) shape
        _anglesOfRotation = np.reshape(_anglesOfRotation, (4,2))
        # Reshape the arrays to the desired (4,2,16) shape. Here, every value in the array is an array of length[16] instead of an float
        _returnAnglesOfRotation = np.reshape(_returnAnglesOfRotation, (4,2,16))

        for _i, _angle in enumerate(_anglesOfRotation):
            # The direction of the direction_axis is the following ;
            # d1 = v_direction1 * directions[0]
            # d2 = v_direction2 * directions[1]

            # Get quaternion for the rotation
            _quaternion1 = MATH.get_quaternion_custom_axis(v_direction1, _angle[0])
            _quaternion2 = MATH.get_quaternion_custom_axis(v_direction2, _angle[1])

            # Rotate the nucleoside by the plane of the base
            _rotatedNuc1 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternion1, nuc1, nuc1[index_origin[0]])
            _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternion2, nuc2, nuc2[index_origin[1]])

            # Calculate the rotation
            _storedDistances[_i] = MATH.get_length_of_vector(_rotatedNuc1[idxDistanceBetweenSubsequentNucs[0]], _rotatedNuc2[idxDistanceBetweenSubsequentNucs[1]])


        _storedDistances = MATH.smallest_difference(_storedDistances, 1.6)
        _indexMin = np.where(_storedDistances == _storedDistances.min())[0][0]
        return _returnAnglesOfRotation[_indexMin][0], _returnAnglesOfRotation[_indexMin][1]
        #v_directions = direction_of_rotation[index_min[0]][0]
        #return v_direction1 * v_directions[0], v_direction2 * v_directions[1]



def assert_rotation_of_bases_by_angle(array_nucs : list, v_directions : list,  index_origin : Union[list, int], _idxAngleBetweenSubsequentNucs: list, angle_to_fit : float)  -> np.ndarray:
    """ Assert the direction of the rotation axis by approximating the best fitting angle between the two nucleic acids; by the backbone"""
    ## ASSERT THE ROTATIONS
    # use a small angle of 2 degrees
#    angle_of_rotation = 5 * (np.pi/180)
    nuc1 = array_nucs[0]
    nuc2 = array_nucs[1]

    idx0 = _idxAngleBetweenSubsequentNucs[0]
    idx1 = _idxAngleBetweenSubsequentNucs[1]
    idx2 = _idxAngleBetweenSubsequentNucs[2]

    # Initialise the values which we use for the ranges of rotation angles over which we rotate over
    _rot0 = np.deg2rad(5)
    _rotI = np.deg2rad(355)
    _rot00, _rot01, _rotI0, _rotI1 = np.deg2rad(2.5), np.deg2rad(20), np.deg2rad(357.5), np.deg2rad(360)

    #angles_of_rotation = np.array([5 * deg2rad, (355 * deg2rad)])
    #v_direction = v_directions[0]
    ##direction_of_rotation = np.array([[1], [-1]])
    #_storedDistances = np.zeros(2)

    # One direction to check
    if len(v_directions) == 1:
        #direction_of_rotation = np.array([[1], [-1]])
        _storedAngles = np.empty(2)
        _anglesOfRotation = np.array([_rot0, _rotI])
        v_direction = v_directions[0]
        _returnAnglesOfRotation = np.array([np.linspace(_rot00, _rot01, 16),
                                           np.linspace(_rotI0, _rotI1, 16)])


        for _i, _angle in enumerate(_anglesOfRotation):
            # The direction of the direction_axis is the following ;
            #d1 = v_direction * directions

            # Get quaternion for the rotation
            _quaternionNuc2 = MATH.get_quaternion_custom_axis(v_direction, _angle)

            # Rotate the nucleoside by the plane of the base
#            _tmpNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, nuc2, nuc2[index_angles])
#            _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, nuc2, nuc2[index_origin])
            _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, nuc2, nuc2[index_origin])

            # Normalise vectors
            v0 = MATH.return_normalized(_rotatedNuc2[idx0] - _rotatedNuc2[idx1])
            v1 = MATH.return_normalized(nuc1[idx2] - _rotatedNuc2[idx1])

            # Calculate the rotation by the angle that the vectors form. The output of the angle is in radians
            _storedAngles[_i] = MATH.get_angle_of_rotation(v0, v1)

        # Find the angle with the lowest difference possible and use that specific set of rotations
        _storedAngles = MATH.smallest_difference(_storedAngles, angle_to_fit)
        _indexMin = np.where(_storedAngles == _storedAngles.min())[0][0]
        #direction = direction_of_rotation[index_min[0]][0]
        return _returnAnglesOfRotation[0] if _indexMin == 0 else _returnAnglesOfRotation[1]

        #return v_direction1 * direction

    # Two sets of directions to check
    if len(v_directions) == 2:
        v_direction1, v_direction2 = v_directions[0], v_directions[1]
        _storedAngles = np.zeros(4)
        #direction_of_rotation = np.array([[1,1], [1,-1], [-1, 1], [-1, -1]])
        _anglesOfRotation = np.array([_rot0, _rot0,
                                      _rot0, _rotI,
                                      _rotI, _rot0,
                                      _rotI, _rotI])

        _returnAnglesOfRotation = np.array([np.linspace(_rot00, _rot01, 16), np.linspace(_rot00, _rot01, 16),
                                            np.linspace(_rot00, _rot01, 16), np.linspace(_rotI0, _rotI1, 16),
                                            np.linspace(_rotI0, _rotI1, 16), np.linspace(_rot00, _rot01, 16),
                                            np.linspace(_rotI0, _rotI1, 16), np.linspace(_rotI0, _rotI1, 16)])

        # Reshape the arrays to the desired (4,2) shape
        _anglesOfRotation = np.reshape(_anglesOfRotation, (4,2))
        # Reshape the arrays to the desired (4,2,16) shape. Here, every value in the array is an array of length[16] instead of an float
        _returnAnglesOfRotation = np.reshape(_returnAnglesOfRotation, (4,2,16))

#        _storedAngles = np.empty(4)
#        v_direction1, v_direction2 = v_directions[0], v_directions[1]
#        direction_of_rotation = np.array([[1,1], [1,-1], [-1, 1], [-1, -1]])

        for _i, _angle in enumerate(_anglesOfRotation):
            # The direction of the direction_axis is the following ;
            #d1 = v_direction1 * directions[0]
            #d2 = v_direction2 * directions[1]

            # Get quaternion for the rotation
            _quaternionNuc1 = MATH.get_quaternion_custom_axis(v_direction1, _angle[0])
            _quaternionNuc2 = MATH.get_quaternion_custom_axis(v_direction2, _angle[1])

            # Rotate the nucleoside by the plane of the base
            _rotatedNuc1 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc1, nuc1, nuc1[index_origin[0]])
            _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, nuc2, nuc2[index_origin[1]])
#            _tmpNuc1 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc1, nuc1, nuc1[index_origin[0]])
#            _tmpNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, nuc2, nuc2[index_origin[1]])

            # Normalise vectors
            v0 = MATH.return_normalized(_rotatedNuc2[idx0] - _rotatedNuc2[idx1])
            v1 = MATH.return_normalized(_rotatedNuc1[idx2] - _rotatedNuc2[idx1])

            # Calculate the rotation by the angle that the vectors form. The output of the angle is in radians
            _storedAngles[_i] = MATH.get_angle_of_rotation(v0, v1)

        # Find the angle with the lowest difference possible and use that specific set of rotations
        _storedAngles = MATH.smallest_difference(_storedAngles, angle_to_fit)
        _indexMin = np.where(_storedAngles == _storedAngles.min())[0][0]
        return _returnAnglesOfRotation[_indexMin][0], _returnAnglesOfRotation[_indexMin][1]
        #v_directions = direction_of_rotation[index_min[0]][0]

        #return v_direction1 * v_directions[0], v_direction2 * v_directions[1]


def distance_or_angle(distances_arr : np.array, angles_arr : np.array, ref_distance : float, ref_angle : float) -> int:
    """ Check to see if the nucleoside fits best distance-bounds or the angle-bounds. Return the index of the array where it is the smallest. """
    stored_distances = MATH.smallest_difference(distances_arr, ref_distance)
    stored_angles = MATH.smallest_difference(angles_arr, ref_angle)

    # Relative differences
    rel_distance = stored_distances.min() / ref_distance
    rel_angle = stored_angles.min() / ref_angle

    # Pick an index from the most suitable one
    if rel_angle < rel_distance:
        return np.where(stored_angles == stored_angles.min())[0]
    return np.where(stored_distances == stored_distances.min())[0]


def tilt_array_to_get_a_better_fit(compl_nuc, compl_linker, prev_nuc, prev_linker, compl_nuc_arr : np.ndarray, complementary_strand : np.ndarray, index_compl : int) -> np.ndarray :
    """ Have the returned array from assert_possible_base_conformations fit the leading strand better. Here we will turn the return array according to length, then dihedral
        compl_nuc is a json object
        prev_nuc is a json object   """

    ## Now we tilt the array slighted if necessary
    # The 'if necessary' will be depending on the length and the dihedral it makes

    # Atom Parsing List for the knowing which atoms to parse from the respective arrays ; for bond length and dihedral evaluation
    APL = PARSE.Atom_Parsing_List(compl_nuc, compl_linker, prev_nuc)
    dihedral = compl_nuc.get_dihedral("alpha")
    angle_to_fit = compl_nuc.get_angle("alpha")

    # Parse the last atom needed from the complementary strand. NB : index_compl has been set to the correct value in the function ' assert_possible ... _and_fit() '.
    # This is the last final atom in the backbone ( 3'-> 5), so in the case of DNA, this would store the vector of the phosphorus atom in its linker
    idx_compl_strand = PARSE.retrieve_atom_index(prev_nuc, APL[3]) + index_compl + prev_linker.mol_length
    v_compl_strand = complementary_strand[idx_compl_strand]

    # Parse the other atoms needed from the current nucleotide that is being fitted
    idx_v0 = PARSE.retrieve_atom_index(compl_nuc, APL[0]) + compl_linker.mol_length
    idx_v1 = PARSE.retrieve_atom_index(compl_nuc, APL[1]) + compl_linker.mol_length
    idx_v2 = PARSE.retrieve_atom_index(compl_linker, APL[2])

    v0 = compl_nuc_arr[idx_v0]
    v1 = compl_nuc_arr[idx_v1]
    v2 = compl_nuc_arr[idx_v2]

    # Calculate of the distance. We will check later if it is between 1 and 2 Angstrom.
    calculated_length = MATH.get_length_of_vector(v2, v_compl_strand)

    # Calculate the dihedral. Let's say the dihedral should not deviate more than 25 degrees?
    calculated_dihr = MATH.dihedral_single(v0, v1, v2, v_compl_strand)

    # Create booleans to assert the calculated values and see if they are within ranges of concordance
    bool_length = MATH.assert_length_of_vector(calculated_length)

    bool_dihedral = MATH.assert_dihedral_of_nucleotide(calculated_dihr, dihedral)

    # If both values are true, just return the array as is, since it is already relatively well positioned
    if bool_length and bool_dihedral :
        return compl_nuc_arr

    # If one or both values are false, try to position the bases to an appropriate orientation
    else :
        ## Generate a quaternion that orients the alpha dihedral properly
        #angle_of_rot_TEST = 5 * (np.pi/180)
        # Retrieve which the denominator of the base that the nucleotide is
        compl_base = compl_nuc.get_base_denominator()

######################################################################################################################### HOW DO WE ASSIGN THE DIRECTION AXIS TO ROTATE ON
        # Make a custom direction and angle. Not to worry, the direction is perpendicular to the movement
        # First try
#        compl_base_atoms1, _ = LCst.retrieve_atoms_for_plane_rotation_of_complement(compl_base, compl_base)
#        compl_base_atoms2, _ , __= LCst.retrieve_atoms_for_positioning_of_complement1(compl_base, compl_base)
#        idx_compl_base_atoms1 = PARSE.retrieve_atom_index(compl_nuc, compl_base_atoms1[0]) + compl_linker.mol_length
#        idx_compl_base_atoms2 = PARSE.retrieve_atom_index(compl_nuc, compl_base_atoms2[2]) + compl_linker.mol_length
#        v_direction = MATH.return_normalized(compl_nuc_arr[idx_compl_base_atoms2] - compl_nuc_arr[idx_compl_base_atoms1])
        # Second try
        complNucleobaseAtoms1, _ , _= CONSTANTS.retrieve_atoms_for_positioning_of_complement1(compl_base, compl_base)
        complNucleobaseAtom2 = CONSTANTS.retrieve_atom_for_direction_axis(compl_base)
        _idxDirNucleobaseAtom1 = PARSE.retrieve_atom_index(compl_nuc, complNucleobaseAtoms1[2]) + compl_linker.mol_length
        _idxDirNucleobaseAtom2 = PARSE.retrieve_atom_index(compl_nuc, complNucleobaseAtom2) + compl_linker.mol_length
        v_direction = MATH.return_normalized(compl_nuc_arr[_idxDirNucleobaseAtom1] - compl_nuc_arr[_idxDirNucleobaseAtom2])

        # angles of rotation that need to be operated with
        #array_of_rotation_angles = np.linspace(2.5, 30, 16) * (np.pi/180)
        #array_of_inverted_rotation_angles = np.linspace(360 - 2.5, 360 - 30, 16) * (np.pi/180)

        ### DISTANCE
        _idxDistanceBetweenSubsequentNucs = [PARSE.retrieve_atom_index(compl_linker, APL[2]),
                                    PARSE.retrieve_atom_index(prev_nuc, APL[3]) + index_compl + prev_linker.mol_length]
#        idx_distance_between_nuc = [PARSE.retrieve_atom_index(compl_linker, APL[2]),
#                                    PARSE.retrieve_atom_index(prev_nuc, APL[3])]

        _anglesOfRotation = assert_rotation_of_bases_by_distance([complementary_strand, compl_nuc_arr], [v_direction], _idxDirNucleobaseAtom1, _idxDistanceBetweenSubsequentNucs)
#        _anglesOfRotation = assert_rotation_of_bases_by_distance([complementary_strand, compl_nuc_arr], [v_direction], _idxDirNucleobaseAtom2, _idxDistanceBetweenSubsequentNucs)
#        v_direction = assert_rotation_of_bases_by_distance([complementary_strand, compl_nuc_arr], [v_direction], _idxDirNucleobaseAtom2, idxDistanceBetweenSubsequentNucs)
#        v_direction = assert_rotation_of_bases_by_distance([complementary_strand, compl_nuc_arr], [v_direction], idx_compl_base_atoms1, idx_distance_between_nuc)
#        v_direction = assert_rotation_of_bases_by_distance([prev_nuc.array, compl_nuc_arr], [v_direction], idx_compl_base_atoms1, idx_distance_between_nuc)

        _storedDistances = np.zeros(shape=(_anglesOfRotation.shape[0]))
        for _i in range(len(_anglesOfRotation)) :
            _quaternionDistance = MATH.get_quaternion_custom_axis(v_direction, _anglesOfRotation[_i])
            testnuc = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionDistance, compl_nuc_arr, compl_nuc_arr[_idxDirNucleobaseAtom1])
#            testnuc = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionDistance, compl_nuc_arr, compl_nuc_arr[_idxDirNucleobaseAtom2])
#            testnuc = MATH.move_to_origin_ROTATE_move_back_to_loc(nucD_quaternion, compl_nuc_arr, compl_nuc_arr[idx_compl_base_atoms1])
            _storedDistances[_i] = MATH.get_length_of_vector(testnuc[idx_v2], complementary_strand[idx_compl_strand])

            if MATH.assert_length_of_vector(_storedDistances[_i]):
#                print(_anglesOfRotation[_i])
                return testnuc

        #print("Accessing angles")
        ### ANGLE
        _idxAngleBetweenSubsequentNucs = [PARSE.retrieve_atom_index(compl_nuc, APL[1]) + compl_linker.mol_length,
                                 PARSE.retrieve_atom_index(compl_linker, APL[2]),
                                 PARSE.retrieve_atom_index(prev_nuc, APL[3]) + index_compl + prev_linker.mol_length]
#        idx_angle_between_nuc = [PARSE.retrieve_atom_index(compl_nuc, APL[1]) + compl_linker.mol_length,
#                                 PARSE.retrieve_atom_index(compl_linker, APL[2]),
#                                 PARSE.retrieve_atom_index(prev_nuc, APL[3])]

        _anglesOfRotation = assert_rotation_of_bases_by_angle([complementary_strand, compl_nuc_arr], [v_direction], _idxDirNucleobaseAtom1, _idxAngleBetweenSubsequentNucs, angle_to_fit)
#        _anglesOfRotation = assert_rotation_of_bases_by_angle([complementary_strand, compl_nuc_arr], [v_direction], _idxDirNucleobaseAtom2, _idxAngleBetweenSubsequentNucs, angle_to_fit)
#        v_direction = assert_rotation_of_bases_by_angle([complementary_strand, compl_nuc_arr], [v_direction], _idxDirNucleobaseAtom2, idxDistanceBetweenSubsequentNucs, angle_to_fit)
#        v_direction = assert_rotation_of_bases_by_angle([complementary_strand, compl_nuc_arr], [v_direction], idx_compl_base_atoms1, idx_angle_between_nuc, angle_to_fit)
#        v_direction = assert_rotation_of_bases_by_angle([prev_nuc.array, compl_nuc_arr], [v_direction], idx_compl_base_atoms1, idx_angle_between_nuc, angle_to_fit)
        idx0 = _idxAngleBetweenSubsequentNucs[0]
        idx1 = _idxAngleBetweenSubsequentNucs[1]
        idx2 = _idxAngleBetweenSubsequentNucs[2]

        _storedAngles = np.zeros(shape=(_anglesOfRotation.shape[0]))
        for _i, _angle in enumerate(_anglesOfRotation) :
            _quaternionAngle = MATH.get_quaternion_custom_axis(v_direction, _angle)
            testnuc = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionAngle, compl_nuc_arr, compl_nuc_arr[_idxDirNucleobaseAtom1])
#            testnuc = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionAngle, compl_nuc_arr, compl_nuc_arr[_idxDirNucleobaseAtom2])
#            testnuc = MATH.move_to_origin_ROTATE_move_back_to_loc(nucA_quaternion, compl_nuc_arr, compl_nuc_arr[idx_compl_base_atoms1])

            # Normalise vectors
            v0 = MATH.return_normalized(testnuc[idx0] - testnuc[idx1])
            v1 = MATH.return_normalized(complementary_strand[idx2] - testnuc[idx1])
#            v1 = MATH.return_normalized(prev_nuc.array[idx2] - testnuc[idx1])

            # Calculate the rotation by the angle that the vectors form. The output of the angle is in radians
            _storedAngles[_i] = MATH.get_angle_of_rotation(v0, v1)

            ## If the angle is suitable according to the boundaries, return the nucleoside's array
            if MATH.assert_size_of_angle(_storedAngles[_i], angle_to_fit):
                return testnuc


        #print("Keep flat")
        # If the angle or distance assertion is not working, just return it as is
        # This returns the nucleotide that has not been fitted; it is just on the plane of the leading strand's nucleobase
        #return compl_nuc_arr


        # Check to see if distance or angle fits the best of our needs
        _indexMin = distance_or_angle(_storedDistances, _storedAngles, 1.6, angle_to_fit)
        _angleOfRot = _anglesOfRotation[_indexMin]
        nuc1_quaternion = MATH.get_quaternion_custom_axis(v_direction, float(_angleOfRot))

        return MATH.move_to_origin_ROTATE_move_back_to_loc(nuc1_quaternion, compl_nuc_arr, compl_nuc_arr[_idxDirNucleobaseAtom2])
#        return MATH.move_to_origin_ROTATE_move_back_to_loc(nuc1_quaternion, compl_nuc_arr, compl_nuc_arr[idx_compl_base_atoms1])


def parse_indexes_of_the_array_for_linker_reorientation(APL : list, nuc, link, nucPN, idx : int = 0) -> Union[np.ndarray, np.ndarray]:
        """ Refactor the part of the code where we parse the correct indexes of the atoms we need from the array. """

        nuc_atom_idx = PARSE.retrieve_atom_index(nuc, APL[1], (idx + link.mol_length + nucPN.mol_length ))
        centrallink_atom_idx = PARSE.retrieve_atom_index(link, APL[2], (idx + nucPN.mol_length))
        nucPN_atom_idx = PARSE.retrieve_atom_index(nucPN, APL[3], idx)
        link_atoms_idx = PARSE.retrieve_atom_index_MULTIPLE(link, link.atom_list, (idx + nucPN.mol_length))

        # Remove central atom from the link_atoms_idx array, since this will make it easier to calculate with later on without making mistakes
        centroid_idx = PARSE.retrieve_atom_index(link, APL[2])
        link_atoms_idx = np.delete(link_atoms_idx, centroid_idx)

        return np.array([nuc_atom_idx, centrallink_atom_idx, nucPN_atom_idx]), link_atoms_idx

#        if ID == "lead":
#            # Parse all the indexes from the AtomParsingList
#            nuc_atom_idx = PARSE.retrieve_atom_index(nuc, APL[1], (idx + link.mol_length + nucPN.mol_length ))
#            centrallink_atom_idx = PARSE.retrieve_atom_index(link, APL[2], (idx + nucPN.mol_length))
#            nucPN_atom_idx = PARSE.retrieve_atom_index(nucPN, APL[3], idx)
#            link_atoms_idx = PARSE.retrieve_atom_index_MULTIPLE(link, link.atom_list, (idx + nucPN.mol_length))
#
#            # Remove central atom from the link_atoms_idx array, since this will make it easier to calculate with later on without making mistakes
#            centroid_idx = PARSE.retrieve_atom_index(link, APL[2])
#            link_atoms_idx = np.delete(link_atoms_idx, centroid_idx)
#
#            return np.array([nuc_atom_idx, centrallink_atom_idx, nucPN_atom_idx]), link_atoms_idx
#
#
#        if ID == "complementary":
#            # Parse all the indexes from the AtomParsingList
#            nuc_atom_idx = PARSE.retrieve_atom_index(nuc, APL[1], (idx + link.mol_length + nucPN.mol_length))
#            centrallink_atom_idx = PARSE.retrieve_atom_index(link, APL[2], (idx + nucPN.mol_length))
#            nucPN_atom_idx= PARSE.retrieve_atom_index(nucPN, APL[3], idx)
#
#            link_atoms_idx = PARSE.retrieve_atom_index_MULTIPLE(link, link.atom_list, (idx + nucPN.mol_length))
#
#            # Remove central atom from the link_atoms_idx array, since this will make it easier to calculate with later on without making mistakes
#            centroid_idx = PARSE.retrieve_atom_index(link, APL[2])
#            link_atoms_idx = np.delete(link_atoms_idx, centroid_idx)
#
#            return np.array([nuc_atom_idx, centrallink_atom_idx, nucPN_atom_idx]), link_atoms_idx



def assert_and_reorient_the_position_of_the_linker(ArrOfIdxBB: np.ndarray, ArrOfIdxLink: np.ndarray, SA : np.ndarray) -> Union[bool, np.ndarray] :
    """ Assess whether or not it is useful to reorient the array of the linker moiety.
        SA is the strand's array

        Important is that the direction of rotation, aka the cross product, is readily normalised 
        To orient the linker properly, we take the vector from the backbone of the nucleotides (O5' and O3') and make the cross product of the plane
            in which the linker lies in parallel to the backbone vector. This ensures that the plane of the linker is perpendicular to that of the vector itself.
        This removes any possible clashes between the linker moiety and the adjacent nucleic acids. """

    # Create the plane in which the backbone atoms lie; e.g. a plane with [O5', P, O3']
    va0 = MATH.return_normalized(SA[ArrOfIdxBB[0]] - SA[ArrOfIdxBB[1]])
    va1 = MATH.return_normalized(SA[ArrOfIdxBB[2]] - SA[ArrOfIdxBB[1]])
    aCross = MATH.get_direction_of_rotation(va1, va0)

    # Create the plane in which the linker's atoms lie; e.g. a plane with [P, OP1, OP2]
    vb0 = MATH.return_normalized(SA[ArrOfIdxLink[0]] - SA[ArrOfIdxBB[1]])
    vb1 = MATH.return_normalized(SA[ArrOfIdxLink[1]] - SA[ArrOfIdxBB[1]])
    bCross = MATH.get_direction_of_rotation(vb1, vb0)

    # An offset of 15 degrees to suggest that over this point, we will start a rotation
    deg2rad = np.pi/180
    NinetyDegrees = 1.5707963267948966                          # 90 degrees in radians
    offSetRadians = 15 * deg2rad
    dotProductInRad = np.arccos(np.dot(aCross, bCross))

    # if the two cross-vectors are more or less orthogonal, return the function
    if NinetyDegrees - offSetRadians <= dotProductInRad <= NinetyDegrees + offSetRadians:
        return SA[ArrOfIdxLink]

    # Generate a cone vector where the angle between aCross and the cone is 90deg.
    cone_vector = MATH.generate_cone_vector(NinetyDegrees)
    quaternion_cone = MATH.get_quaternion(aCross * -1.0)
    rotated_cone_vector = MATH.rotate_with_quaternion(quaternion_cone, cone_vector)
    # The vectors in the cone will be evaluated against the vector that goes from O5' -> O3' 
    parallelVector = MATH.return_normalized(SA[ArrOfIdxBB[2]] - SA[ArrOfIdxBB[0]])
    quaternion = MATH.get_quaternion(parallelVector, bCross)

    #return MATH.move_to_origin_ROTATE_move_back_to_loc(quaternion, SA[ArrOfIdxLink], SA[ArrOfIdxBB[1]])
    return MATH.move_to_origin_ROTATE_move_back_to_loc(quaternion, SA[ArrOfIdxLink], SA[ArrOfIdxBB[1]])



def capping_retrieve_atomnames(list_of_leading_sequence : list, list_of_complementary_sequence : list, backbone_dict : dict, nuc_dict : dict) -> list:
    """ This should return the names of the capping atoms.
            i.e. ; if the H binds to O5', the atom name of the H should be HO5' """
    # Initiate a list with filler values
    list_lead_atom = [0, 0]
    list_compl_atom = [0, 0]
    index_list = [-1, 0]
    for atom in range(len(index_list)):
        # Leading strand
        lead_nucleoside_fname = nuc_dict[list_of_leading_sequence[index_list[atom]]][0]
        atom_lead = initMolecule.Nucleoside(lead_nucleoside_fname)
        # Complementary strand
        compl_nucleoside_fname = nuc_dict[list_of_complementary_sequence[index_list[atom]]][0]
        atom_compl = initMolecule.Nucleoside(compl_nucleoside_fname)

        backbone_lead = backbone_dict[atom_lead.get_nucleic_acid_code()]
        backbone_compl = backbone_dict[atom_compl.get_nucleic_acid_code()]

        # Parse the atom of interest from the backbone
        list_lead_atom[atom] = backbone_lead[index_list[atom]]
        list_compl_atom[atom] = backbone_compl[index_list[atom]]

    # Concatenate the list
    list_of_atoms = list_lead_atom + list_compl_atom

    # Correct the atom names
    for i in range(len(list_of_atoms)):
        list_of_atoms[i] = ["H" + list_of_atoms[i]]

    return list_of_atoms


def capping_retrieve_atomarrays(leading_array : np.ndarray, list_of_leading_sequence : list, complementary_array : np.ndarray, list_of_complementary_sequence : list, backbone_dict : dict, nuc_dict : dict) -> np.ndarray:
    """ This should return the coordinates of the capping atoms.
        We are going to make a set angle and dihedral for all the hydrogens, since it's not going to matter that much in which angle they are.
            As long as there are no clashes and the dihedral and angle are within ranges of reason. """

    # Hardcode the angle and dihedral, since they just have to be out of the way and not cause clashes.
    angle = 135 * (np.pi/180)
    dihedral = 179

    # Initialise the list of atoms to parse from
    lead_atom_list = np.zeros(shape=(2,3), dtype=object)
    compl_atom_list = np.zeros(shape=(2,3), dtype=object)

    # Initialise the array of the supposed vector
    H_vectors = np.zeros(shape=(4,3), dtype=object)

    # Initialise a dictionary that parses the atomnames by index
    index_list = [0, -1]

    # Parse the correct atoms from the backbone dict. This for loop prints the keys of the dict, when calling the variable 'atom'
    for i in range(len(index_list)):
        # Leading strand
        lead_nucleoside_fname = nuc_dict[list_of_leading_sequence[index_list[i]]][0]
        lead_atom = initMolecule.Nucleoside(lead_nucleoside_fname)
        # Complementary strand
        compl_nucleoside_fname = nuc_dict[list_of_complementary_sequence[index_list[i]]][0]
        compl_atom = initMolecule.Nucleoside(compl_nucleoside_fname)

        backbone_lead = backbone_dict[lead_atom.get_nucleic_acid_code()]
        backbone_compl = backbone_dict[compl_atom.get_nucleic_acid_code()]

        # Now retrieve the index of the arrays for which we want to parse the atoms
        if index_list[i] == 0 :
            # Parse the atom of interest from the backbone
            for j in [-1, -2, -3]:
                lead_atom_list[i][j] = backbone_lead[j]
                compl_atom_list[i][j] = backbone_compl[j]

            lead_atom_indexes = PARSE.retrieve_atom_index_MULTIPLE(lead_atom, lead_atom_list[i])
            compl_atom_indexes = PARSE.retrieve_atom_index_MULTIPLE(compl_atom, compl_atom_list[i])

            # Now retrieve the coordinates from the array
            l0 = leading_array[lead_atom_indexes[-3]]
            l1 = leading_array[lead_atom_indexes[-2]]
            l2 = leading_array[lead_atom_indexes[-1]]

            c0 = complementary_array[compl_atom_indexes[-3]]
            c1 = complementary_array[compl_atom_indexes[-2]]
            c2 = complementary_array[compl_atom_indexes[-1]]

            # Now calculate the would be dihedral and retrieve the H position in xyz
            vector_lead = UB.generate_vector_of_interest(angle, dihedral, [l2, l1, l0])
            vector_compl = UB.generate_vector_of_interest(angle, dihedral, [c2, c1, c0])

            # Correct the length, add the vector to the array and store it in H_vectors
            H_vectors[0] = (MATH.return_normalized(vector_lead) * 1.2) + l2
            H_vectors[2] = (MATH.return_normalized(vector_compl) * 1.2) + c2

        elif index_list[i] == -1 :
            for j in [0, 1, 2]:
                lead_atom_list[i][j] = backbone_lead[j]
                compl_atom_list[i][j] = backbone_compl[j]

            # Create the index counter for both lead and compl strand
            lead_index = leading_array.shape[0] - lead_atom.mol_length
            compl_index = complementary_array.shape[0] - compl_atom.mol_length

            lead_atom_indexes = PARSE.retrieve_atom_index_MULTIPLE(lead_atom, lead_atom_list[i], lead_index)
            compl_atom_indexes = PARSE.retrieve_atom_index_MULTIPLE(compl_atom, compl_atom_list[i], compl_index)

            # Now retrieve the coordinates from the array
            l0 = leading_array[lead_atom_indexes[0]]
            l1 = leading_array[lead_atom_indexes[1]]
            l2 = leading_array[lead_atom_indexes[2]]

            c0 = complementary_array[compl_atom_indexes[0]]
            c1 = complementary_array[compl_atom_indexes[1]]
            c2 = complementary_array[compl_atom_indexes[2]]

            # Now calculate the would be dihedral and retrieve the H position in xyz
            vector_lead = UB.generate_vector_of_interest(angle, dihedral, [l0, l1, l2])
            vector_compl = UB.generate_vector_of_interest(angle, dihedral, [c0, c1, c2])

            # Correct the length, add the vector to the array and store it in H_vectors
            H_vectors[1] = (MATH.return_normalized(vector_lead) * 1.2) + l0
            H_vectors[3] = (MATH.return_normalized(vector_compl) * 1.2) + c0


    return H_vectors


