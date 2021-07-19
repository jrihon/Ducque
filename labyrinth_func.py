import numpy as np
import pandas as pd
import json, sys
from typing import Union

import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2
import labyrinth_func_tools3 as LFT3


"""
Labyrinth_func.py
The script that contains the classes and all the functions that concatenate the workflow of consecutively adding the linker and nucleotides.
"""

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
        return json.loads(self.jason["identity"])[3][0]

class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_shape(self) -> int:
        """ returns the size of the linker shape, but only the first value """
        return int(json.loads(self.jason["pdb_properties"]["Shape"])[0])




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
    APL = ["C4'", "C5'", "O5'", "P", "OP2", "OP1"]

    # Retrieve the vectors of the atoms that make up the dihedral you research
    # Dihedral C4' - C5' - O5' - P
    id_v0 = LFT2.retrieve_atom_index(nucleoside, APL[0])
    v0 = nucl_array[id_v0]

    id_v1 = LFT2.retrieve_atom_index(nucleoside, APL[1])
    v1 = nucl_array[id_v1]

    id_v2 = LFT2.retrieve_atom_index(nucleoside, APL[2])
    v2 = nucl_array[id_v2]

    # Find the vector that corresponds to the O5' -> P vector
    single_vector1 = generate_vector_of_interest(nucleoside.get_angle("beta"), nucleoside.get_dihedral("beta"), [v2, v1, v0])

    # Add the single_vector1 to O5' atom, to denote the location of P. Then move the linker array to that location
    link = LFT1.position_phosphate(v2, single_vector1, linker.array)   # the linker is in the correct position  

    ## We will rotate the linker twice, to get the correct orientation of the linker in 3D space
    # Dihedral C5' - O5' - P - OP2
    id_v3 = LFT2.retrieve_atom_index(linker, APL[3])
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
    leading_strand : the nucleic acid strand to which we append next_nucleoside to
    """
    # The first thing to do is to find the location of the subsequent atom, here O3', then rotate the next_nucleoside by zeta and epsilon
    # Afterwards, we turn over the epsilon dihedral and by rotate the normal of the plane have to the plane we want; this positions everything!

    # Atom Parsing List (ATP) = Parse which linker and which nucleotide the previous one is
    APL = ["C5'", "O5'", "P", "O3'", "C3'", "C4'"]
    # Dihedral Parsing List (DPL) = Parse which dihedrals are required to rotate on and over
    DPL = ["alpha", "zeta", "epsilon"]
    # Angle Parsing List (AngPL)
    AngPL = ["alpha", "zeta", "epsilon"]
    # Get the indices of the vectors you do not want moved
    exclusion_list = ["O3'", "C3'"]

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
            next_nucleoside_originloc_N = LFT1.move_vector_to_origin(next_nucleoside_loc, distance_to_origin_N)
            # exclude atoms from the rotation and perform the rotation of the planes
            exclude_these_atoms = [LFT2.retrieve_atom_index(next_nucleoside, atom) for atom in exclusion_list]
            next_nucleoside_originloc_rotated = LFT1.apply_rotation_of_planes(quaternion_plane, next_nucleoside_originloc_N, exclude_these_atoms)
            next_nucleoside_loc = LFT1.move_vector_to_loc(next_nucleoside_originloc_rotated, distance_to_origin_N)
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

        next_nucleoside_originloc_N = LFT1.move_vector_to_origin(next_nucleoside_loc, distance_to_origin_N)
        next_nucleoside_originloc_rotated = LFT1.rotate_with_quaternion(quaternion_N, next_nucleoside_originloc_N)
        next_nucleoside_loc = LFT1.move_vector_to_loc(next_nucleoside_originloc_rotated, distance_to_origin_N)


def position_complementary_base(leading_base, complementary_base, leading_array : np.ndarray, index_counter : int) -> np.ndarray:
    """ This functions positions the base correctly, which inherently positions the entire nucleoside correctly.
        After this function, there comes an iterative process of fitting the backbone correctly.

        leading_base is a json object
        complementary_base is a json object                     """

    ### DATA PARSING
    ## Parse the correct atoms to get started for
    # Parse the type of base from the nucleoside json object
    base1 = leading_base.get_base_denominator()
    base2 = complementary_base.get_base_denominator()

    rotPlane_base1_atoms, rotPlane_base2_atoms = LFT3.retrieve_atoms_for_plane_rotation_of_complement(base1, base2)

    ## Get the proper vectors for the plane rotation
    # leading base vectors
    list_rotPlane_base1_id = LFT2.retrieve_atom_index_MULTIPLE(leading_base, rotPlane_base1_atoms, index_counter)

    vL0 = leading_array[list_rotPlane_base1_id[0]]
    vL1 = leading_array[list_rotPlane_base1_id[1]]
    vL2 = leading_array[list_rotPlane_base1_id[2]]

    pL0 = LFT1.return_normalized(vL1 - vL0)     # C2 - N1
    pL1 = LFT1.return_normalized(vL2 - vL0)     # C6 - N1

    # complementary base vectors
    list_rotPlane_base2_id = LFT2.retrieve_atom_index_MULTIPLE(complementary_base, rotPlane_base2_atoms)

    vC0 = complementary_base.array[list_rotPlane_base2_id[0]]
    vC1 = complementary_base.array[list_rotPlane_base2_id[1]]
    vC2 = complementary_base.array[list_rotPlane_base2_id[2]]

    pC0 = LFT1.return_normalized(vC1 - vC0)     # C4 - N9
    pC1 = LFT1.return_normalized(vC2 - vC0)     # C8 - N9

### Rotate the plane
    ## Rotate the complementary base onto the plane of the leading base
    # Get the cross product of the atoms of leading_base and reverse the sign of that vector. Dont forget to normalise
    #cross_base1 = LFT1.return_normalized(np.cross(pL0, pL1))
    cross_base1 = LFT1.return_normalized(np.cross(pL1, pL0))
    #cross_base2 = LFT1.return_normalized(np.cross(pC0, pC1))
    cross_base2 = LFT1.return_normalized(np.cross(pC1, pC0))

    # Get the cross product of the atoms of leading_base that make it if both bases are in the same plane, the cross product vectors are exactly opposite
    plane1_quaternion = LFT1.get_quaternion( cross_base1 * -1.0, cross_base2)

    # Apply the rotation
    base2_at_origin = LFT1.move_vector_to_origin(complementary_base.array, vC0)
    base2_at_origin_rotated = LFT1.rotate_with_quaternion(plane1_quaternion, base2_at_origin)

### Apply a translation to get the initial positioning
    ### (X1) Now we have rotated the plane, we calculate X1 vector and move our complementary base there.
    compl1_base1_atoms, compl1_base2_atom = LFT3.retrieve_atoms_for_positioning_of_complement1(base1, base2)

    list_compl1_base1_id = LFT2.retrieve_atom_index_MULTIPLE(leading_base, compl1_base1_atoms, index_counter)

    # Get the required parameters for the rotations
    Q_angle, Q_dihedral, R_angle, R_dihedral = LFT3.retrieve_angles_and_dihedrals_for_initial_base_positioning(base1)
    dist_between_bases1 = 1.81 ; dist_between_bases2 = 1.87

    # Calculate for X1 then and position Hn (compl_base) to the position
    single_vector_compl1 = generate_vector_of_interest(Q_angle, Q_dihedral, [leading_array[list_compl1_base1_id[2]], leading_array[list_compl1_base1_id[1]], leading_array[list_compl1_base1_id[0]]])
    single_vector_compl1 = LFT1.return_normalized(single_vector_compl1) * dist_between_bases1

    ## Now that we have the Hn position, we move the complementary base to the position by doing the following :
    # Get distance at which we want to position the base's hydrogen bond
    move_to_this_loc1 = leading_array[list_compl1_base1_id[2]] + single_vector_compl1
    # retrieve the vector of the atom that we want to attach 'move_to_this_loc' to.
    v_compl1_base2_id = LFT2.retrieve_atom_index(complementary_base, compl1_base2_atom)

    # the distance from the location we want to have compl1_base2_atom move to. But we move the array from the origin, so substract that and then move it
    distance_for_compl1 = move_to_this_loc1 - base2_at_origin_rotated[v_compl1_base2_id]

    base2_at_loc = LFT1.move_vector_to_loc(base2_at_origin_rotated, distance_for_compl1)

### Apply rotation a second time
    ### (X2) Now we have placed complementary base in the correct orientation, now we rotate it to get the 3D orientation in order by turning onto the second hydrogen bond
    compl2_base1_atoms, compl2_base2_atom = LFT3.retrieve_atoms_for_position_of_complement2(base1, base2)
    list_compl2_base1_id = LFT2.retrieve_atom_index_MULTIPLE(leading_base, compl2_base1_atoms, index_counter)

    # Calculate for X2
    single_vector_compl2 = generate_vector_of_interest(R_angle, R_dihedral, [leading_array[list_compl2_base1_id[2]], leading_array[list_compl2_base1_id[1]], leading_array[list_compl2_base1_id[0]]])
    single_vector_compl2 = LFT1.return_normalized(single_vector_compl2) * dist_between_bases2

    # Creates the coordinate at which we want to move our compl2_base2_atom to
    move_to_this_loc2 = leading_array[list_compl2_base1_id[2]] + single_vector_compl2

    v_compl2_base2_id = LFT2.retrieve_atom_index(complementary_base, compl2_base2_atom)
    # j1 = MVTL2 - Hn
    # j2 = Ot - Hn
    # move the 
    v_Hn = base2_at_loc[v_compl1_base2_id]
    j1 = LFT1.return_normalized(move_to_this_loc2 - v_Hn)
    j2 = LFT1.return_normalized(base2_at_loc[v_compl2_base2_id] - v_Hn)

    plane2_quaternion = LFT1.get_quaternion(j1, j2)

    base2_at_origin = base2_at_loc - base2_at_loc[v_compl1_base2_id]
    base2_at_origin_rotated = LFT1.rotate_with_quaternion(plane2_quaternion, base2_at_origin)

    base2_at_loc = base2_at_origin_rotated + v_Hn

    return base2_at_loc





def generate_complementary_sequence(sequence_list : list, complement : Union[list, str]) -> list:
    """ sequence list is the given input.
        complement will specify what the complementary strand will look like. choices between homo - DNA - RNA """
    complementary_dictDNA = { "A" : "T", "T" : "A", "G" : "C", "C" : "G", "U" : "A" }
    complementary_dictRNA = { "A" : "U", "T" : "A", "G" : "C", "C" : "G", "U" : "A" }

    bases = LFT2.retrieve_bases_list(sequence_list)
    if isinstance(complement, list):
        complementary_sequence = list(map(lambda x: x.strip(","), complement))

        # See of the lengths match. If the lengths do not match, give an assertion error and print the following string. This exits the software too.
        assert len(sequence_list) == len(complementary_sequence), "The length of the complementary strand does not match the length of the leading strand!"

        return complementary_sequence

    if complement.lower() == "homo":
        chemistry = LFT2.retrieve_chemistry_list(sequence_list)

        # Switch the bases the get their complementary base
        if chemistry[0] == "r":
            comp_bases = LFT2.get_complementary_bases(bases, complementary_dictRNA)
            complementary_sequence = LFT2.concatenate_chem_and_bases(chemistry, comp_bases)
            return complementary_sequence

        # If the complementary strand is not RNA based, use the DNA native bases.
        comp_bases = LFT2.get_complementary_bases(bases, complementary_dictDNA)
        complementary_sequence = LFT2.concatenate_chem_and_bases(chemistry, comp_bases)
        return complementary_sequence

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


def create_PDB_from_array_final(leading_array : np.ndarray, list_of_leading_sequence : list, complementary_array : np.ndarray, list_of_complementary_sequence : list) -> None:
    """ Write out the data for the pdb filename
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html"""

    list_of_leading_sequence = list_of_leading_sequence[::-1]
    #list_of_complementary_sequence = list_of_complementary_sequence[::-1]

    # LEADING STRAND
    df_leading = pd.DataFrame()

    df_leading["RecName"] = ["ATOM" for x in range(leading_array.shape[0])]
    df_leading["AtomNum"] = np.arange(start=1, stop=leading_array.shape[0] + 1)
    df_leading["AtomName"] = LFT2.LEAD_pdb_AtomNames_or_ElementSymbol(list_of_leading_sequence, "Atoms")
    df_leading["AltLoc"] = " "
    df_leading["ResName"] = LFT2.LEAD_pdb_Residuename(list_of_leading_sequence)
    df_leading["Chain"] = "A"
    df_leading["Sequence"] = LFT2.LEAD_pdb_Sequence(list_of_leading_sequence)
    df_leading["X_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,0]))
    df_leading["Y_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,1]))
    df_leading["Z_coord"] = list(map(lambda x: "{:.3f}".format(x), leading_array[:,2]))
    df_leading["Occupancy"] = "1.00"
    df_leading["Temp"] = "0.00"
    df_leading["SegmentID"] = str("   ")
    df_leading["ElementSymbol"] = LFT2.LEAD_pdb_AtomNames_or_ElementSymbol(list_of_leading_sequence, "Symbol")

    # TER line between the single strands
    TER_line = pd.DataFrame({"RecName" : "TER", "AtomNum" : 69}, index=[69])
    ln_lead = len(df_leading["AtomNum"])

    # COMPLEMENTARY STRAND
    df_complementary = pd.DataFrame()

    df_complementary["RecName"] = ["ATOM" for x in range(complementary_array.shape[0])]
    df_complementary["AtomNum"] = np.arange(start=ln_lead + 1, stop=(ln_lead + complementary_array.shape[0] + 1))
    df_complementary["AtomName"] = LFT2.COMPLEMENTARY_pdb_AtomNames_or_ElementSymbol(list_of_complementary_sequence, "Atoms")
    df_complementary["AltLoc"] = " "
    df_complementary["ResName"] = LFT2.COMPLEMENTARY_pdb_Residuename(list_of_complementary_sequence)
    df_complementary["Chain"] = "K"
    df_complementary["Sequence"] = LFT2.COMPLEMENTARY_pdb_Sequence(list_of_complementary_sequence, len(list_of_leading_sequence))
    df_complementary["X_coord"] = list(map(lambda x: "{:.3f}".format(x), complementary_array[:,0]))
    df_complementary["Y_coord"] = list(map(lambda x: "{:.3f}".format(x), complementary_array[:,1]))
    df_complementary["Z_coord"] = list(map(lambda x: "{:.3f}".format(x), complementary_array[:,2]))
    df_complementary["Occupancy"] = "1.00"
    df_complementary["Temp"] = "0.00"
    df_complementary["SegmentID"] = str("   ")
    df_complementary["ElementSymbol"] = LFT2.COMPLEMENTARY_pdb_AtomNames_or_ElementSymbol(list_of_complementary_sequence, "Symbol")

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
                pdb.write("%-6s%5s%5s%s%2s%3s%5d  %8s%8s%9s%6s%7s%4s     %2s\n" % tuple(split_line))
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
