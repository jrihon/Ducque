import sys, os, json
import numpy as np

import labyrinth_func as LabF
import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2

""" Create dictionary of the filename and their molecule code """

codex_acidum_nucleicum = {
'dT': ['json/dna_thymidine.json', 'json/phosphate_linker.json'],
'dA': ['json/dna_adenosine.json', 'json/phosphate_linker.json'],
'dC': ['json/dna_cytosine.json', 'json/phosphate_linker.json'],
'dG': ['json/dna_guanosine.json', 'json/phosphate_linker.json'],


}

def printexit(*arg):
    """Function to print out any numbers of parameters and then exit the software
    This is for debugging purposes	"""
    print(arg)

    sys.exit(0)


def Architecture(nucleic_acid_list):

    ## PARSE DATA FROM JSON
    nucleic_acid = nucleic_acid_list[0]

    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
    nucleo = LabF.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

    # Parse dictionary for the correct linker segment
    linker = LabF.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

    # Retrieve the vectors of the atoms that make up the dihedral you research
    id_O5 = LFT2.retrieve_atom_index(nucleo, "O5'")
    v2 = nucleo.array[id_O5]

    id_C5 = LFT2.retrieve_atom_index(nucleo, "C5'")
    v1 = nucleo.array[id_C5]

    id_C4 = LFT2.retrieve_atom_index(nucleo, "C4'")
    v0 = nucleo.array[id_C4]

    single_vector1 = LabF.generate_vector_of_interest(linker.get_COP(), nucleo.get_beta(), [v2, v1, v0])

    # POSITION THE LINKER TO THE LOCATION THE NEW VECTOR POINTS AT
    # Adding single vector, which points to P to 
    link = LFT1.position_linker(v2, single_vector1, linker)

#------------ SO FAR THE LINKER HAS BEEN POSITIONED, BUT IT NEEDS TO ROTATE  
    # Dihedral C5' - O5' - P - OP2
    id_P = LFT2.retrieve_atom_index(linker, "P")
    v3 = link[id_P]

    id_OP2 = LFT2.retrieve_atom_index(linker, "OP2")
    v4 = link[id_OP2]

    single_vector2 = LabF.generate_vector_of_interest(linker.get_OPO2(), linker.get_OP2_dihedral(), [v3, v2, v1])

    # The distance from P to the origin of the xyz system
    p_to_origin = v3

    # move the linker to the origin, by positioning the posphorus at [0,0,0]
    link_to_origin = link - p_to_origin
#------------ Rotate the phosphate group
    ## Rotate the phosphate for a first time
    # Define the vector that goes from P to OP2 and normalize it
    P_OP2 = LFT1.return_normalized(v4 - v3)

    # get quaternion to rotate the linker a first time
    quaternion_P1 = LFT1.get_quaternion(single_vector2 , vector_to_rotate_from=P_OP2)
    #Rotate it
    link = LFT1.rotate_with_quaternion(quaternion_P1, link_to_origin)
    # Move the rotated linker to the place of origin
    link = link + p_to_origin
#------------ ROTATE THE LINKER AGAIN, BUT ON THE AXIS OF P_OP2, to get P_OP1
    # Dihedral C5' - O5' - P - OP1
    # Since we have rotated P -> OP2, we need to override the vector again from the array 'link' we just overrided
    v4 = link[id_OP2]

    id_OP1 = LFT2.retrieve_atom_index(linker, "OP1")
    v5 = link[id_OP1]

    single_vector3 = LabF.generate_vector_of_interest(linker.get_OPO1(), linker.get_OP1_dihedral(), [v3, v2, v1])

#------------ Rotate OP_1 of the linker.array onto it single_vector3
    # OP2 - P ; P -> OP2
    P_O2 = LFT1.return_normalized(v4 - v3)
    # OP1 - P ; P -> OP1
    P_O1 = LFT1.return_normalized(v5 - v3)

    # generate quaternion
    quaternion_P2 = LFT1.get_quaternion_custom_axis(single_vector3, P_O1, P_O2 * -1.0)
    # move linker back to origin
    link = link - p_to_origin
    # rotate the vector
    link = LFT1.rotate_with_quaternion(quaternion_P2, link)
    # Bring OP_final to the location of the phosphorus
    link = link + p_to_origin

######################### LINKER POSITIONED ###############################################

######################### POSITION THE NEXT NUCLEOTIDE ####################################
# Start with a the stacked array, that makes it easier
    nucleotide = np.vstack((link, nucleo.array))

    # import the next nucleoside
    nextnuc_acid = nucleic_acid_list[1]

    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
    nextnuc = LabF.Nucleoside(codex_acidum_nucleicum[nextnuc_acid][0])

#### We position the nextnuc by the position of O3' and then rotate it correctly by the zeta dihedral

#------------ Position nextnuc correctly
    # dihedral C5' (v1) - O5' (v2) - P(v3) - O3'. All vectors are already present, so no make new variables

    single_vector4 = LabF.generate_vector_of_interest(linker.get_O3PO5(), nextnuc.get_alpha(), [v3, v2, v1])

    # Add single_vector4 to , so that we define the location of O3'
    v6 = LFT1.move_vector_to_loc(LFT1.return_normalized(single_vector4) * 1.6, v3)
    # O3' index
    id_O3 = LFT2.retrieve_atom_index(nextnuc, "O3'")
    # Get distance from nextnuc O3' to the position defined as O3'
    O3_distance_from_origin = nextnuc.array[id_O3]
    nextnuc_distance = v6 - O3_distance_from_origin
    # Get the distance from the P_O3 now and add it to nextnuc_origin
    #   which will move the entire molecule(nextnuc) to the position of O3'
    nextnuc_loc = LFT1.move_vector_to_loc(nextnuc.array, nextnuc_distance)

#------------ Rotate nextnuc correctly along zeta
    # O5' (v2), P (v3), O3' (v6), C3'
    P_O3_C3 = 119.032 * (np.pi / 180)
    zeta_dihr = nucleo.get_zeta()

    single_vector5 = LabF.generate_vector_of_interest(P_O3_C3, zeta_dihr, [v6, v3, v2])

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
    epsilon_dihr = nucleo.get_epsilon()
    # Get angle O3' - C3' - C4'
    O3_C3_C4 = 111.919 * (np.pi/180)

    single_vector6 = LabF.generate_vector_of_interest(O3_C3_C4, epsilon_dihr, [v7, v6, v3])

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


######################### CREATE THE PDB THAT GOES WITH ARRAY INPUTTED ####################
    stacked_array = (nextnuc_loc, link, nucleo.array)
    LabF.create_PDB_from_matrix(np.vstack(stacked_array), nucleic_acid_list)


########################                    FINAL                    ###################### 
