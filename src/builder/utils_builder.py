import json, sys
from typing import Union, List

import initMolecule

import numpy as np
import builder.mathematics as MATH
import builder.parse_or_write as PARSE
import builder.builder_constants as CONSTANTS
import builder.utils_of_utils_builder as utilsUB
from ducquelib.library import TABLE_NUCLEOTIDES, TABLE_CONFORMATIONS, TABLE_BACKBONE

import systemsDucque as SD

TN = TABLE_NUCLEOTIDES
TC = TABLE_CONFORMATIONS
TB = TABLE_BACKBONE


""" utils_builder.py
The script that contains the classes and all the functions that concatenate the workflow of consecutively adding the linker and nucleotides. """


def generate_vector_of_interest(angle : float, dihedral : float, atom_array : list) -> np.ndarray:
    """ This function generates a single vector for which there exists only one angle and dihedral.
        The atom_array contains the three first atoms in the sequence that make up the dihedral.
        Example = if the sequence of a dihedral is C4' - C5' - O5' - P (beta backbone), then the atom_array is [O5', C5', C4'] 
        atom_array : list[np.ndarray]                                                                                           """

    # Get the vector to rotate the cone vector onto. This is done on the middle two atoms of the sequence
    # Example : [O5'] minus (-) [C5'] results in a vector that goes [C5' -> O5']
    p0 = MATH.return_normalized(atom_array[0] - atom_array[1])

    # Generate the cone vector
    cone_vector = MATH.generate_cone_vector(angle)

    # Get the quaternion that corresponds to the desired rotation
    # Note, if 'vector_to_from" is empty, it defaults to the Z-axis, around which the cone is generated np.array([0,0,1])
    quaternion = MATH.get_quaternion(p0 * -1.0)

    # Rotate the cone vector
    rotated_cone_vector = MATH.rotate_with_quaternion(quaternion, cone_vector)

    # Calculate the possible dihedrals we have with the first three atoms and the cone vector. Remember: the cone vector represents the position of the fourth atom
    range_of_dihedrals = MATH.dihedral_array(atom_array, rotated_cone_vector)

    # Now extrapolate in between the range_of_dihedrals to get the correct theta angle to generate the single vector of interest.
    # The theta angle corresponds to the angle (spherical coordinates r, theta, phi) at which we need to generate our vector at.
    theta_interpolate = MATH.get_interpolated_dihedral(range_of_dihedrals, dihedral)

    # Generate the correct vector of interest
    single_vector = MATH.generate_and_rotate_single_vector(theta_interpolate, angle, quaternion)

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
    APL = PARSE.Atom_Parsing_List(nucleoside, linker)

    # Retrieve the vectors of the atoms that make up the dihedral you research
    # Dihedral C4' - C5' - O5' - P
    id_v0 = PARSE.retrieve_atom_index(nucleoside, APL[0])
    id_v1 = PARSE.retrieve_atom_index(nucleoside, APL[1])
    id_v2 = PARSE.retrieve_atom_index(nucleoside, APL[2])
    id_v3 = PARSE.retrieve_atom_index(linker, APL[3])

    v0 = nucl_array[id_v0]
    v1 = nucl_array[id_v1]
    v2 = nucl_array[id_v2]

    # Find the vector that corresponds to the O5' -> P vector
    single_vector1 = generate_vector_of_interest(nucleoside.get_angle("beta"), nucleoside.get_dihedral("beta"), [v2, v1, v0])

    # Position the phosphate by moving the linker's center atom to the origin and then bringing to the location that was calculated for
    move_link_to = v2 + (single_vector1 * 1.6)
    link_array = MATH.move_vector_to_loc(linker.array - linker.array[id_v3], move_link_to)

    ## We will rotate the linker twice, to get the correct orientation of the linker in 3D space
    # Dihedral C5' - O5' - P - OP2
    v3 = link_array[id_v3]

    single_vector2 = generate_vector_of_interest(linker.get_angle("angle_1"), linker.get_dihedral("dihedral_1"), [v3, v2, v1])

    # This is the distance from the linker's atom to the origin
    link_distance = link_array[id_v3]

    # move the linker to the origin, by positioning the phosphorus at [0,0,0]
    link_at_origin = MATH.move_vector_to_origin(link_array, link_distance)

    # Rotate the linker a first time
    # Define the vector that goes from P to OP2 and normalize it
    id_v4 = PARSE.retrieve_atom_index(linker, APL[4])
    v4 = link_array[id_v4]
    p3_4 = MATH.return_normalized(v4 - v3)

    # Get quaternion to rotate the linker a first time and rotate it
    quaternion_P1 = MATH.get_quaternion(single_vector2 , p3_4)
    rotated_link_at_origin = MATH.rotate_with_quaternion(quaternion_P1, link_at_origin)

    # Move the rotated linker back to the calculated position of the phosphorus atom
    link_array = MATH.move_vector_to_loc(rotated_link_at_origin, link_distance)


    # Dihedral C5' - O5' - P - OP1
    id_v5 = PARSE.retrieve_atom_index(linker, APL[5])
    v5 = link_array[id_v5]

    v4 = link_array[id_v4]
    p3_4 = MATH.return_normalized(v4 - v3)

    # Generate vector we want to rotate P_OP1 on to
    single_vector3 = generate_vector_of_interest(linker.get_angle("angle_1"), linker.get_dihedral("dihedral_2"), [v3, v2, v1])
    p3_5 = MATH.return_normalized(v5 - v3)

    plane_rotate_onto = MATH.get_direction_of_rotation(single_vector2, single_vector3)
    plane_rotate_from = MATH.get_direction_of_rotation(single_vector2, p3_5)

    quaternion_P2 = MATH.get_quaternion(plane_rotate_onto , plane_rotate_from)

    # move the linker to the origin, by positioning the phosphorus at [0,0,0]
    link_at_origin2 = MATH.move_vector_to_origin(link_array, link_distance)

#    # Get quaternion to rotate the linker a first time and rotate it
#    # Rotate only part that needs to be moved
#    quaternion_P2 = MATH.get_quaternion(single_vector3 , p3_5)
#    part_to_rotate = link_at_origin2[id_v5:]
    rotated_link_at_origin2 = MATH.rotate_with_quaternion(quaternion_P2, link_at_origin2)
#
#    link_array_at_origin = np.vstack((link_at_origin2[:id_v5], rotated_part_at_origin2))
#
#    # Move the rotated linker back to the calculated position of the phosphorus atom
    link_array = MATH.move_vector_to_loc(rotated_link_at_origin2, link_distance)
#
    # Stack the arrays on top of each other
    nucleotide = np.vstack((link_array, nucl_array))

    return nucleotide


def generate_complementary_sequence(sequence_list : list, complement : Union[List[str], str]) -> list:
    """ sequence list is the given input.
        complement will specify what the complementary strand will look like. choices between homo - DNA - RNA """

    import process_CLI_inputs
    complementary_dictDNA = { "A" : "T", "T" : "A", "G" : "C", "C" : "G", "U" : "A" }
    complementary_dictRNA = { "A" : "U", "T" : "A", "G" : "C", "C" : "G", "U" : "A" }

    # The keys, meaning the nucleosides, from the complementary dictionary wil be parsed as a list
    keys_of_dict = list(TN.keys())

    # Get a list of the bases of the leading strand, to later build a complementary strand
    bases = PARSE.retrieve_bases_list(sequence_list)

    ## IF THE INPUT OF THE COMPLEMENT IS A LIST
    if isinstance(complement, list):

        complementary_sequence = list(map(lambda x: x.strip(","), complement))

        # See of the lengths match. If the lengths do not match, give an assertion error and print the following string. This exits the software too.
        assert len(sequence_list) == len(complementary_sequence), "The length of the complementary strand does not match the length of the leading strand!"

        # Check if one of the nucleosides in the prompted list is wrong, i.e. not existing or wrongly prompted (misspelled)
        # returns boolean
        if not process_CLI_inputs.check_if_nucleotides_are_valid(keys_of_dict):
            sys.exit(0)

        # Reverse the prompted list
        return complementary_sequence[::-1]


    # IF THE SEQUENCE NEEDS TO BE A HOMODUPLEX
    if complement.upper() == "HOMO":
        complementary_sequence = []
        list_of_chemistries = PARSE.retrieve_chemistry_list(sequence_list)

        # Get only the nucleosides of the same chemistry as the chemistry it will complement
        for chem_i in range(len(list_of_chemistries)):
            # Get the possible nucleosides we can choose from with that have the same chemistry as in the leadstrand's nucleoside
            parsed_chemistry, ln_str = PARSE.retrieve_chemistry(list_of_chemistries[chem_i])
            homo_chemistry_list = PARSE.retrieve_homo_nucleosides(keys_of_dict, parsed_chemistry, ln_str)

            # Assess which base fits the complementary nucleoside
            complementary_nucleoside = PARSE.assess_possible_complementary_base(parsed_chemistry, bases[chem_i], homo_chemistry_list, complementary_dictRNA, complementary_dictDNA)
            complementary_sequence.append(complementary_nucleoside)

        return complementary_sequence


    # IF THE SEQUENCE NEEDS TO BE A HETERODUPLEX; CHECK THE VALIDITY OF THE PROMPTED CHEMISTRY AND RETURN A LIST WITH POSSIBLE ABBREVIATED CHEMISTRY CODES
    NUC_ID = process_CLI_inputs.check_if_chemistry_is_valid(complement)

    # Return the abbreviated name of the chemistry
    chemCode = PARSE.return_chemistrycode(NUC_ID)

    # If the 'chemCode + T' does not exist as a nucleotide, this means that the the uracil variant is valid as the nucleobase of the chemistry
    if chemCode + "T" in keys_of_dict :
        comp_bases = PARSE.get_complementary_bases(bases, complementary_dictDNA)
        complementary_sequence = PARSE.concatenate_chem_and_bases(chemCode, comp_bases)
        return complementary_sequence

    comp_bases = PARSE.get_complementary_bases(bases, complementary_dictRNA)
    complementary_sequence = PARSE.concatenate_chem_and_bases(chemCode, comp_bases)
    return complementary_sequence


def assert_leading_strand_nucleotide_conformation(conformations, prev_nucleoside, prev_link, leading_strand):
    """ Position the different conformations and assert which conformation is the most parralel with respect to the bases of both nucleosides """

    # If there is only one conformation possible, no point in asserting the best one
    if len(conformations) == 1:
        return conformations[0]

    # Get the json object parse the correct atoms for the vector comparison
    next_nucleoside = initMolecule.Nucleoside(conformations[0])

    prev_base = prev_nucleoside.get_base_denominator()
    next_base = next_nucleoside.get_base_denominator()
    prev_base_atoms, next_base_atoms = CONSTANTS.retrieve_atoms_for_plane_rotation_of_complement(prev_base, next_base)

    # Parse the indexes of the atoms of interest
    prev_base_atoms_indexes = PARSE.retrieve_atom_index_MULTIPLE(prev_nucleoside, prev_base_atoms, index_counter=prev_link.mol_length)

    # Previous nucleoside Cross Vector
    cross_prev = utilsUB.return_cross_vector_for_plane_rotation(leading_strand, prev_base_atoms_indexes)

    # Array of different cross vectors of the next nucleoside
    possible_cross_vectors = np.zeros(len(conformations), dtype=object)

    for conf in range(len(conformations)):
        # Previous nucleoside Cross Vector
        nucleoside = initMolecule.Nucleoside(conformations[conf])
        # Position the nucleoside in place
        arrays_of_the_conformations = utilsUB.position_next_nucleoside(nucleoside, prev_nucleoside, prev_link, leading_strand)
        # Parse the indexes of the atoms of interest
        next_base_atoms_indexes = PARSE.retrieve_atom_index_MULTIPLE(next_nucleoside, next_base_atoms)
        possible_cross_vectors[conf] = utilsUB.return_cross_vector_for_plane_rotation(arrays_of_the_conformations, next_base_atoms_indexes)

    # Calculate the scalar product of the two vectors, convert it to radians. The smallest value in the array corresponds to the dot(cross1, cross2) with the least amount of difference
    # In simple terms, it means both planes of the respective bases are parallel
    radians_product_array = np.zeros(len(conformations), dtype=object)

    for v in range(len(conformations)):
        scalar_product = np.dot(cross_prev, possible_cross_vectors[v])
        radians_product_array[v] = np.arccos(scalar_product)

#    stored_distances = MATH.smallest_difference(radians_product_array, 1)
    smallest_distance = radians_product_array.min()
    index_smallest_distance = np.where(radians_product_array == smallest_distance)[0][0]

    return conformations[index_smallest_distance]


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
    lead1_base, _ = initMolecule.Nucleoside(lead_bases[0][0]), initMolecule.Linker(lead_bases[0][1])
    lead2_base, lead2_link = initMolecule.Nucleoside(lead_bases[1][0]), initMolecule.Linker(lead_bases[1][1])

    # Prepare variables for the forloop
    index_lead -= lead1_base.mol_length
    size_of_lead_base2 = lead2_base.mol_length + lead2_link.mol_length

    # This works for for when there is a single conformation or multiple conformations
    for i in range(len(compl1_base_confs)):
        # Instantiate the complementary base
        compl1_base = initMolecule.Nucleoside(compl1_base_confs[i])
        # position the complementary base
        compl1_base_arr = utilsUB.position_complementary_base(lead1_base, compl1_base, leading_strand, index_lead)

        for j in range(len(compl2_base_confs)):
            #decrement the index_lead
            index_lead -= size_of_lead_base2

            # Instantiate the complementary base
            compl2_base = initMolecule.Nucleoside(compl2_base_confs[j])
            # position the complementary base
            compl2_base_arr = utilsUB.position_complementary_base(lead2_base, compl2_base, leading_strand, index_lead)
            # add linker to the complementary base
            compl2_base_arr = position_phosphate_linker(compl2_base, compl2_base_arr, compl2_linker)

            # Store the array of the dinucleotide
            array_of_possible_dinucleotide_conformations[i][j] = np.vstack((compl1_base_arr, compl2_base_arr))

            #increment the index_lead back
            index_lead += size_of_lead_base2


    # decrement it again, since we return this variable for the upcoming nucleotides after this function finalises.
    index_lead -= size_of_lead_base2

    # Parse the index of the values we want. compl1 is the atom that is being attached to the linker of the previous nucleotide. compl2 is the last atom in the backbone of the linker.
    APL1 = PARSE.Atom_Parsing_List(compl2_base, compl2_linker, compl1_base)

    _idxDistCompl1 = PARSE.retrieve_atom_index(compl1_base, APL1[3])
    _idxDistCompl2 = PARSE.retrieve_atom_index(compl2_linker, APL1[2]) + compl1_base.mol_length # Add object.mol_length because we parse from a single array that is a nucleoside

    # Now we iterate over the array. When we calculate for a small distance, we will remember the index of its array and keep that as the best fit so far
    _distance = 100                          # start with an unreasonable distance. This value will be tested against the distances we find
    _idxCompl1, _idxCompl2 = -1, -1       # start with defaulted values. If they do not change, kill the process
    for i in range(len(compl1_base_confs)):
        for j in range(len(compl2_base_confs)):
            tmp_array = array_of_possible_dinucleotide_conformations[i][j]
            v_compl1 = tmp_array[_idxDistCompl1]
            v_compl2 = tmp_array[_idxDistCompl2]

            distance_between_nucleotides = MATH.get_length_of_vector(v_compl1, v_compl2)
            #print(distance_between_nucleotides)
            if 1.00 <= distance_between_nucleotides < _distance :
                # Override the variable dist and save the best dinucleotide conformation as an array. Also save off the index to parse the correct json file.
                _distance = distance_between_nucleotides
                arr_best_fit = array_of_possible_dinucleotide_conformations[i][j]

                _idxCompl1 = i
                _idxCompl2 = j

    # If these variables have not been assigned a value, it probably means the input dihedrals are not optimal
    if not utilsUB.assess_complX_id(_idxCompl1, _idxCompl2):
        sys.exit("Try out a different dihedral torsion for your inputs.\n"
                "With the given input dihedrals/angles, the nucleotides of the complementary strand are very close together\n"
                "This might have to do the shape of the leading strand. ")


    # If the distance is already of the desired size, then return the dinucleotide as is
    if MATH.assert_length_of_vector(_distance):
        complementary_strand = arr_best_fit
        nuc1 = initMolecule.Nucleoside(compl1_base_confs[_idxCompl1])
        nuc2 = initMolecule.Nucleoside(compl2_base_confs[_idxCompl2])
        return complementary_strand, index_lead, nuc1, nuc2


    # Now that we have the best fit, just to be sure we instantiate the correct object for the correct json filem for atom parsing.
    # This is generally not necessary, if the only differences between the json files, of the conformations of the same nucleosides, is the atoms coordinates.
    compl1_nuc = initMolecule.Nucleoside(compl1_base_confs[_idxCompl1])
    compl2_nuc = initMolecule.Nucleoside(compl2_base_confs[_idxCompl2])
    # get the size of the array as an index to correctly parse the required atoms of compl2_base
    index_fit = compl1_nuc.mol_length

    # Parse the two nucleotides from the array
    _arrNuc1 = arr_best_fit[:index_fit]
    _arrNuc2 = arr_best_fit[index_fit:]

    _idxDistCompl1 = PARSE.retrieve_atom_index(compl1_nuc, APL1[3])
    _idxDistCompl2 = PARSE.retrieve_atom_index(compl2_nuc, APL1[1])
    #id_dist_compl2 = PARSE.retrieve_atom_index(compl2_linker, APL1[2])

    ## GET THE VECTORS AND THE INDEX OF THE ATOMS REQUIRED FOR THE ROTATIONS
    # NUCLEOTIDE 1
    _nucNucleobase1 = compl1_base.get_base_denominator()
#    nuc1_origin, _ = CONSTANTS.retrieve_atoms_for_plane_rotation_of_complement(nuc1_base, nuc1_base)
    _nucAtomOrigin1 = CONSTANTS.retrieve_atom_for_direction_axis(_nucNucleobase1)
    _nucAtomsForDirection1, _, _ = CONSTANTS.retrieve_atoms_for_positioning_of_complement1(_nucNucleobase1, _nucNucleobase1)
#    id_nuc1_origin = PARSE.retrieve_atom_index(compl1_nuc, nuc1_origin[0])
    _idxNucAtomsOrigin1 = PARSE.retrieve_atom_index(compl1_nuc, _nucAtomOrigin1)
    _idxNucAtomsForDirection1 = PARSE.retrieve_atom_index(compl1_nuc, _nucAtomsForDirection1[2])

    # NUCLEOTIDE 2
    _nucNucleobase2 = compl2_base.get_base_denominator()
#    nuc2_origin, _ = CONSTANTS.retrieve_atoms_for_plane_rotation_of_complement(nuc2_base, nuc2_base)
    _nucAtomOrigin2 = CONSTANTS.retrieve_atom_for_direction_axis(_nucNucleobase2)
    _nucAtomsForDirection2, _ , _ = CONSTANTS.retrieve_atoms_for_positioning_of_complement1(_nucNucleobase2, _nucNucleobase2)
#    id_nuc2_origin = PARSE.retrieve_atom_index(compl2_nuc, nuc2_origin[0]) + compl2_linker.mol_length
    _idxNucAtomsOrigin2 = PARSE.retrieve_atom_index(compl2_nuc, _nucAtomOrigin2) + compl2_linker.mol_length
    _idxNucAtomsForDirection2 = PARSE.retrieve_atom_index(compl2_nuc, _nucAtomsForDirection2[2]) + compl2_linker.mol_length

    # Initialise direction axis
    v_direction1 = MATH.return_normalized(_arrNuc1[_idxNucAtomsForDirection1] - _arrNuc1[_idxNucAtomsOrigin1])
    v_direction2 = MATH.return_normalized(_arrNuc2[_idxNucAtomsForDirection2] - _arrNuc2[_idxNucAtomsOrigin2])

    # Make lists of the variables, just because it makes it easier to read
    array_nucs = [_arrNuc1, _arrNuc2]
#    index_origin = [_idxNucAtomsOrigin1, _idxNucAtomsOrigin2]                                      ########## CHECK
    index_origin = [_idxNucAtomsForDirection1, _idxNucAtomsForDirection2]                           ########## CHECK
    index_distance_between_nuc = [_idxDistCompl1, _idxDistCompl2]
    v_directions = [v_direction1, v_direction2]

    # Generate angles of rotation. Make note that we won't go past a 15 degree angle, since that general will not be necessary to turn that much.
    #array_of_rot_angles = np.linspace(2.5, 15, 16) * (np.pi/180)

    ### CALCULATIONS BASED ON DISTANCE
    _arrayOfRotationAngles1, _arrayOfRotationAngles2 = utilsUB.assert_rotation_of_bases_by_distance(array_nucs, v_directions, index_origin, index_distance_between_nuc)

    _storedDistances = np.empty(shape=len(_arrayOfRotationAngles1))
    for i, _ in enumerate(_arrayOfRotationAngles1) :
        _quaternionNuc1 = MATH.get_quaternion_custom_axis(v_direction1, _arrayOfRotationAngles1[i])
        _quaternionNuc2 = MATH.get_quaternion_custom_axis(v_direction2, _arrayOfRotationAngles2[i])

        _rotatedNuc1 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc1, _arrNuc1, _arrNuc1[_idxNucAtomsForDirection1])
        _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, _arrNuc2, _arrNuc2[_idxNucAtomsForDirection2])

        _storedDistances[i] = MATH.get_length_of_vector(_rotatedNuc1[_idxDistCompl1], _rotatedNuc2[_idxDistCompl2])
        # If the distance is suitable according to the boundaries, return the stack
        if MATH.assert_length_of_vector(_storedDistances[i]) :
            complementary_strand = np.vstack((_rotatedNuc1, _rotatedNuc2))
            return complementary_strand, index_lead, compl1_nuc, compl2_nuc



    ### CALCULATIONS BASED ON ANGLE
    idx_angle_between_nuc = [PARSE.retrieve_atom_index(compl2_nuc, APL1[1]) + compl2_linker.mol_length,
                             PARSE.retrieve_atom_index(compl2_linker, APL1[2]),
                             PARSE.retrieve_atom_index(compl1_nuc, APL1[3])]
    angle_to_fit = compl2_base.get_angle("alpha")

    _arrayOfRotationAngles1, _arrayOfRotationAngles2 = utilsUB.assert_rotation_of_bases_by_angle(array_nucs, v_directions, index_origin, idx_angle_between_nuc, angle_to_fit)

    idx0 = idx_angle_between_nuc[0]
    idx1 = idx_angle_between_nuc[1]
    idx2 = idx_angle_between_nuc[2]

    _storedAngles = np.empty(shape=len(_arrayOfRotationAngles2))
    for _i, _ in enumerate(_arrayOfRotationAngles2) :
        _quaternionNuc1 = MATH.get_quaternion_custom_axis(v_direction1, _arrayOfRotationAngles1[_i])
        _quaternionNuc2 = MATH.get_quaternion_custom_axis(v_direction2, _arrayOfRotationAngles2[_i])

        _rotatedNuc1 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc1, _arrNuc1, _arrNuc1[_idxNucAtomsForDirection1])
        _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, _arrNuc2, _arrNuc2[_idxNucAtomsForDirection2])

        # Normalise vectors
        v0 = MATH.return_normalized(_rotatedNuc2[idx0] - _rotatedNuc2[idx1])
        v1 = MATH.return_normalized(_rotatedNuc1[idx2] - _rotatedNuc2[idx1])

        # Calculate the rotation by the angle that the vectors form. The output of the angle is in radians
        _storedAngles[_i] = MATH.get_angle_of_rotation(v0, v1)

        if MATH.assert_size_of_angle(_storedAngles[_i], angle_to_fit) :
            complementary_strand = np.vstack((_rotatedNuc1, _rotatedNuc2))
            return complementary_strand, index_lead, compl1_nuc, compl2_nuc

    # If there is no solution, take the one with the smallest angle
    _storedDistances = MATH.smallest_difference(_storedDistances, angle_to_fit)
    # If this part of the code is reached, that means none of the evaluated distances are suitable. Let's search the shortest one in the list
    _indexMin = np.where(_storedDistances == _storedDistances.min())[0]
    _angleOfRot = _arrayOfRotationAngles2[_indexMin]
    _quaternionNuc1 = MATH.get_quaternion_custom_axis(v_direction1, float(_angleOfRot))
    _quaternionNuc2 = MATH.get_quaternion_custom_axis(v_direction2, float(_angleOfRot))

#    nuc1 = MATH.move_to_origin_ROTATE_move_back_to_loc(quat_nuc1, arr_nuc1, arr_nuc1[id_nuc1_origin])
#    nuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(quat_nuc2, arr_nuc2, arr_nuc2[id_nuc2_origin])
    _rotatedNuc1 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc1, _arrNuc1, _arrNuc1[_idxNucAtomsOrigin1])
    _rotatedNuc2 = MATH.move_to_origin_ROTATE_move_back_to_loc(_quaternionNuc2, _arrNuc2, _arrNuc2[_idxNucAtomsOrigin2])
    complementary_strand = np.vstack((_rotatedNuc1, _rotatedNuc2))
    return complementary_strand, index_lead, compl1_nuc, compl2_nuc


def assert_possible_base_conformations_and_fit(leading_nuc, leading_array : np.ndarray, conformations : list, compl_linker, complementary_strand : np.ndarray,
                                                                                prev_compl_nuc, prev_compl_linker, index_lead : int, index_compl : int) : # -> Nucleoside
    """ Conformations is a list of the TC[NA], that contains one, two or three different conformations of the same NA.

        Conformations contains the names of particular json files, which will be converted to json objects.
        leading_nuc is already a json object.
        compl_linker is already a json object
        prev_compl_nuc is already json objects

        after the initial fit is done, it will calculate the best orientation for the nucleoside to be in, with respect to leading and complementary strand.
        Returns the conformation that has the best fit into the leading strand and the progressing complementary strand """

    # If there is only one conformation, there is no need to assert the differences in distance or dihedral
    if len(conformations) == 1:
        conf_n = initMolecule.Nucleoside(conformations[0])
        conf_n_arr = utilsUB.position_complementary_base(leading_nuc, conf_n, leading_array, index_lead)
        conf_with_linker = position_phosphate_linker(conf_n, conf_n_arr, compl_linker)
        conf_n.array = utilsUB.tilt_array_to_get_a_better_fit(conf_n, compl_linker, prev_compl_nuc, prev_compl_linker, conf_with_linker, complementary_strand, index_compl)
        return conf_n

    ## Since there are multiple conformations available, we will need to sort out which one will fit the best.
    ## First we parse the correct indexes for the vectors we want to evaluate (check the distance)
    ## We first check if there are nucleotides that fit well without rotation of the base-plane. If that's not the case, we will check all conformations based on best fit after rotation of the base-plane.
    ## We then assert the different conformations after they've been fitted to the best it can

    # nuc_data is used just to be able to parse the correct indexes in the code below
    nuc_data = initMolecule.Nucleoside(conformations[0])

    # The first atom in the backbone of the nucleoside the current nucleoside attaches to
    atomOfInterest2 = TB[json.loads(prev_compl_nuc.jsonObject["identity"])[1]][0]   # first value of the TB to parse. else make list(_VAL)[0]
    bb_id = PARSE.retrieve_atom_index(prev_compl_nuc, atomOfInterest2) + index_compl + prev_compl_linker.mol_length
    bb_v = complementary_strand[bb_id]

    # The last atom in the backbone of the linker of the current nucleoside. 
    atomOfInterest1 = TB[json.loads(compl_linker.jsonObject["identity"])[0]][-1] # since it is the last of the TB that should be parsed. else make list(_VAL)[-1]
    link_id = PARSE.retrieve_atom_index(compl_linker, atomOfInterest1)

    ## First we assert the conformations without tilting by the base-plane

    # Instantiate a matrix to store possible conformations in an array
    possibilities_of_conformations = np.zeros(shape=len(conformations), dtype=object)
    # Instantiate a matrix to store the distances in the backbone, from linker to previous nucleotide, and whether or not they fit within the boundaries of asserting set distances
    stored_bb_distances = np.zeros(shape=len(conformations), dtype=object)
    stored_bb_distance_bools = np.zeros(shape=len(conformations), dtype=bool)
    stored_bb_dihedral_bools = np.zeros(shape=len(conformations), dtype=bool)


    for file_n in range(len(possibilities_of_conformations)):
        conf_n = initMolecule.Nucleoside(conformations[file_n])
        conf_n_arr = utilsUB.position_complementary_base(leading_nuc, conf_n, leading_array, index_lead)
        # add the linker to the nucleoside
        possibilities_of_conformations[file_n] = position_phosphate_linker(conf_n, conf_n_arr, compl_linker)

        link_v = possibilities_of_conformations[file_n][link_id]
        # Check the distance between the two nucleotides, so P -> O3' (with native nucs as example) one of which will have the shortest distance. Pick that one
        stored_bb_distances[file_n] = MATH.get_length_of_vector(link_v, bb_v)
        stored_bb_distance_bools[file_n] = MATH.assert_length_of_vector(stored_bb_distances[file_n])
        stored_bb_dihedral_bools[file_n] = utilsUB.assert_the_dihedral_of_interest(conf_n, possibilities_of_conformations[file_n], compl_linker, prev_compl_nuc, complementary_strand, index_compl, prev_compl_linker)

    # check dihedral suitability
    if np.any(stored_bb_dihedral_bools == True) :
        index_of_best_conformation = utilsUB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_dihedral_bools)
        if not stored_bb_distance_bools[index_of_best_conformation] == True:
            conf_n = initMolecule.Nucleoside(conformations[index_of_best_conformation])
            conf_n.array = utilsUB.tilt_array_to_get_a_better_fit(conf_n, compl_linker, prev_compl_nuc, prev_compl_linker, possibilities_of_conformations[index_of_best_conformation], complementary_strand, index_compl)
            return conf_n

    # if the check dihedral suitability fails, check for distance suitability
    if np.any(stored_bb_distance_bools == True):
        index_of_best_conformation = utilsUB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_distance_bools)
        conf_n =  initMolecule.Nucleoside(conformations[index_of_best_conformation])
        conf_n.array = possibilities_of_conformations[index_of_best_conformation]
        return conf_n


    ## Alas, none of the conformations were suitable enough, it seems we will have to fit them over and over.
    # Instantiate a matrix to store possible conformations in an array
    possibilities_of_conformations = np.zeros(shape=len(conformations), dtype=object)
    # Instantiate a matrix to store the distances in the backbone, from linker to previous nucleotide, and whether or not they fit within the boundaries of asserting set distances
    stored_bb_distances = np.zeros(shape=len(conformations), dtype=object)
    stored_bb_distance_bools = np.zeros(shape=len(conformations), dtype=bool)
    stored_bb_dihedral_bools = np.zeros(shape=len(conformations), dtype=bool)

    for file_n in range(len(possibilities_of_conformations)):
        conf_n = initMolecule.Nucleoside(conformations[file_n])
        conf_n_arr = utilsUB.position_complementary_base(leading_nuc, conf_n, leading_array, index_lead)
        # add the linker to the nucleoside
        conf_n_arr_link = position_phosphate_linker(conf_n, conf_n_arr, compl_linker)
        # tilt the conformation to the best possible fit
        possibilities_of_conformations[file_n] = utilsUB.tilt_array_to_get_a_better_fit(conf_n, compl_linker, prev_compl_nuc, prev_compl_linker, conf_n_arr_link, complementary_strand, index_compl)

        link_v = possibilities_of_conformations[file_n][link_id]
        # Check the distance between the two nucleotides, so P -> O3' (with native nucs as example) one of which will have the shortest distance. Pick that one
        stored_bb_distances[file_n] = MATH.get_length_of_vector(link_v, bb_v)
        stored_bb_distance_bools[file_n] = MATH.assert_length_of_vector(stored_bb_distances[file_n])
        stored_bb_dihedral_bools[file_n] = utilsUB.assert_the_dihedral_of_interest(conf_n, possibilities_of_conformations[file_n], compl_linker, prev_compl_nuc, complementary_strand, index_compl, prev_compl_linker)

    if np.any(stored_bb_dihedral_bools == True) :
        index_of_best_conformation = utilsUB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_dihedral_bools)
        conf_n =  initMolecule.Nucleoside(conformations[index_of_best_conformation])
        conf_n.array = possibilities_of_conformations[index_of_best_conformation]
        return conf_n

    # This gives the one with the least distance from the desired bb_distance
    index_of_best_conformation = utilsUB.retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_distance_bools)
    conf_n =  initMolecule.Nucleoside(conformations[index_of_best_conformation])
    conf_n.array = possibilities_of_conformations[index_of_best_conformation]
    return conf_n


def orient_the_linker_moieties_better(CONF_LIST : list, LINK_LIST : list, leading_array : np.ndarray, compl_array : np.ndarray) -> Union[np.ndarray, np.ndarray]:
    """ There are, for N nucleic acids in a single strands, always N - 1 nucleotides and therefor N - 1 linker moieties built in each strand.
        This function is implemented to consider the placement of all the linker moieties. If they are unfortunately placed, moved them to a more suitable conformations.

        Important to note is that we build the leading strand from the bottom up.
        Important to note is that we build the complementary strand from the top up.                         


        EDIT : This function is at the moment not in use.

        """

    # LEADING STRAND
#    lead_nuc = CONF_LIST[0]
#    lead_link = LINK_LIST[0]
#
#    idx_lead = leading_array.shape[0]          # Initiate the index counter for the leading strand. 
#    for nucleotide in range(len(lead_nuc) - 1):
#        # Initialise the molecules as json objects
#        nuc = initMolecule.Nucleoside(lead_nuc[nucleotide])
#        link = initMolecule.Linker(lead_link[nucleotide])
#        nextnuc = initMolecule.Nucleoside(lead_nuc[nucleotide + 1])
#
#        # Decrement the idx_lead to be able to parse the atoms from the leading strands's array
#        idx_lead -= (nuc.mol_length + link.mol_length + nextnuc.mol_length)
#
#        # Parse the required atoms from AtomParsingList
#        lead_APL = PARSE.Atom_Parsing_List(nuc, link, nextnuc)
#        ArrOfIdxBB, ArrOfIdxLink = utilsUB.parse_indexes_of_the_array_for_linker_reorientation(lead_APL, nextnuc, link, nuc, idx_lead)
#
#        # Based on the positioning, the linker moiety has or has not been rotated
#        # It is a bit overkill to override the same value (if the reorientation was not necessary), but overriding two indexes at a time is not such a large consumption of time to worry over
#        linker_array = utilsUB.assert_and_reorient_the_position_of_the_linker(ArrOfIdxBB, ArrOfIdxLink, leading_array)
#        leading_array[ArrOfIdxLink] = linker_array
#
#        # Increment for the next cycle
#        idx_lead += nextnuc.mol_length
#
#
#    # COMPLEMENTARY STRAND
#    compl_nuc = CONF_LIST[1]
#    compl_link = LINK_LIST[1]
#
#    idx_link = 0
#    idx_compl = 0                           # Initiate the index counter for the complementary strand.
#
#    for nucleotide in range(1, len(compl_nuc)):
#
#        # Initialise the molecules as json objects
#        nuc = initMolecule.Nucleoside(compl_nuc[nucleotide])
#        link = initMolecule.Linker(compl_link[idx_link])
#        prevnuc = initMolecule.Nucleoside(compl_nuc[nucleotide - 1])
#
#        # Parse the required atoms from AtomParsingList
#        compl_APL = PARSE.Atom_Parsing_List(prevnuc, link, nuc)
#        ArrOfIdxBB, ArrOfIdxLink = utilsUB.parse_indexes_of_the_array_for_linker_reorientation(compl_APL, nuc, link, prevnuc, idx_compl)
#
#        # Based on the positioning, the linker moiety has or has not been rotated
#        # It is a bit overkill to override the same value (if the reorientation was not necessary), but overriding two indexes at a time is not such a large consumption of time to worry over
#        linker_array = utilsUB.assert_and_reorient_the_position_of_the_linker(ArrOfIdxBB, ArrOfIdxLink, compl_array)
#        compl_array[ArrOfIdxLink] = linker_array
#
#        # Decrement the idx_lead to be able to parse the atoms from the leading strands's array
#        idx_compl += (prevnuc.mol_length + link.mol_length)
#
#        idx_link += 1

    return leading_array, compl_array


def cap_nucleic_acid_strands(leading_array : np.ndarray, leading_sequence : list, complementary_array : np.ndarray, complementary_sequence : list) -> Union[np.ndarray, list]:
    """ Cap the nucleic acid strands with a hydrogen, to finish the build of the duplex
        We create two functions in utilsUB that parse both the correct coordinates of the capping atoms and their names """

    # Retrieve the backbone atoms, retrieve the filenames too from the TB and TN tables

    # Defines and returns the cartesian position of the hydrogens 
    atom_array = utilsUB.capping_retrieve_atomarrays(leading_array, leading_sequence, complementary_array, complementary_sequence, TB, TN)

    # Returns the name of the hydrogens that have been defined with the previous function
    atom_names = utilsUB.capping_retrieve_atomnames(leading_sequence, complementary_sequence, TB, TN)

    return atom_array, atom_names


def create_PDB_from_array_final(outfile : str, leadingArray : np.ndarray, listOfLeadingSequence : list, complementaryArray : np.ndarray, listOfComplementarySequence : list):
    """ Write out the PDB formatted file

        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html"""

    # Reverse the sequence of this to parse the correct atoms along the array
    listOfLeadingSequence = listOfLeadingSequence[::-1]

    # Capping of the duplexes
    # First two elements belong to the leading strand, the last two elements belong to the complementary strand. A total of four elements per variable
    terminalArray, terminalAtomNames = cap_nucleic_acid_strands(leadingArray, listOfLeadingSequence, complementaryArray, listOfComplementarySequence)

    # Adjust leading strand's array for the capping
    # For readability, I had added the completion of the array here and not in a backend function
    leadingArray = np.vstack((terminalArray[0], leadingArray, terminalArray[1]))


    # LEADING STRAND
    leadingStrand = initMolecule.PdbInstance()

    leadingStrand.SetAtomNumber(leadingArray)
    leadingStrand.SetAtomName(terminalAtomNames, listOfLeadingSequence, "lead")
    leadingStrand.SetResidueName(listOfLeadingSequence)
    leadingStrand.SetChainLetter("J")
    leadingStrand.SetSequenceNumber(listOfLeadingSequence, "lead")
    leadingStrand.SetAtomArray(leadingArray, "lead")
    leadingStrand.SetElementSymbols(listOfLeadingSequence)


    # Adjust leading strand's array for the capping.
    # For readability, I had added the completion of the array here and not in a backend function
    complementaryArray = np.vstack((terminalArray[2], complementaryArray, terminalArray[3]))

    # COMPLEMENTARY STRAND
    complementaryStrand = initMolecule.PdbInstance()

    complementaryStrand.SetAtomNumber(complementaryArray, len(leadingStrand.AtomNumber))
    complementaryStrand.SetAtomName(terminalAtomNames, listOfComplementarySequence, "complementary")
    complementaryStrand.SetResidueName(listOfComplementarySequence)
    complementaryStrand.SetChainLetter("R")
    complementaryStrand.SetSequenceNumber(listOfComplementarySequence, "complementary")
    complementaryStrand.SetAtomArray(complementaryArray, "complementary")
    complementaryStrand.SetElementSymbols(listOfComplementarySequence)


    # Check if pdb is already appended to the string name
    if not outfile.split('.')[-1] == 'pdb' :
        outfile = outfile + ".pdb"
    
    # Write out the *.pdb file
    with open(outfile, "w") as PDB:
        # Write out Leading Strand
        for idxL in range(len(leadingStrand.AtomNumber)):
            line = [leadingStrand.RecordName,
                    leadingStrand.AtomNumber[idxL],
                    leadingStrand.AtomName[idxL],
                    leadingStrand.ResidueName[idxL],
                    leadingStrand.Chain,
                    leadingStrand.SequenceNumber[idxL],
                    leadingStrand.x[idxL], leadingStrand.y[idxL],leadingStrand.z[idxL],
                    leadingStrand.Occupancy, leadingStrand.TempFactor,
                    leadingStrand.ElementSymbol[idxL] ]

            PDB.write("%-4s  %5d %4s %3s %s%4d    %8s%8s%8s%6s%6s          %2s\n" % tuple(line))

        # Add TER line in between the two strands
        PDB.write("TER\n")

        # Write out Complementary Strand
        for idxC in range(len(complementaryStrand.AtomNumber)):
            line = [complementaryStrand.RecordName,
                    complementaryStrand.AtomNumber[idxC],
                    complementaryStrand.AtomName[idxC],
                    complementaryStrand.ResidueName[idxC],
                    complementaryStrand.Chain,
                    complementaryStrand.SequenceNumber[idxC],
                    complementaryStrand.x[idxC], complementaryStrand.y[idxC],complementaryStrand.z[idxC],
                    complementaryStrand.Occupancy, complementaryStrand.TempFactor,
                    complementaryStrand.ElementSymbol[idxC] ]

            PDB.write("%-4s  %5d %-4s %3s %s%4d    %8s%8s%8s%6s%6s          %2s\n" % tuple(line))

        # Finalize the *.pdb file
        PDB.write("END")

    SD.print_writing(outfile)
