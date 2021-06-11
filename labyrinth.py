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
    nucleoside = LabF.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

    # Parse dictionary for the correct linker segment
    linker = LabF.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

    leading_strand = LabF.position_phosphate_linker(nucleoside, linker)

    # LINKER IS POSITIONED

    # POSITION THE NEXT NUCLEOTIDE
    # import the next nucleoside
    nextnuc_acid = nucleic_acid_list[1]
    nextnuc = LabF.Nucleoside(codex_acidum_nucleicum[nextnuc_acid][0])

    previous_nucleic_acid = nucleic_acid_list[1-1]
    prevnuc, prevlink = LabF.Nucleoside(codex_acidum_nucleicum[previous_nucleic_acid][0]), LabF.Desmos(codex_acidum_nucleicum[previous_nucleic_acid][1])

    nextnuc_positioned = LabF.position_next_nucleotide(nextnuc, prevnuc, prevlink, leading_strand)
#    # Retrieve the vectors of the atoms that make up the dihedral you research
#    id_O5 = LFT2.retrieve_atom_index(nucleoside, "O5'")
#    d2 = nucleoside.array[id_O5]
#
#    id_C5 = LFT2.retrieve_atom_index(nucleoside, "C5'")
#    d1 = nucleoside.array[id_C5]
#
#    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
#    nextnuc = LabF.Nucleoside(codex_acidum_nucleicum[nextnuc_acid][0])
#
##### We position the nextnuc by the position of O3' and then rotate it correctly by the zeta dihedral
#
##------------ Position nextnuc correctly
#    # dihedral C5' (v1) - O5' (v2) - P(v3) - O3'. All vectors are already present, so no make new variables
#    v_P = nucleotide[LFT2.retrieve_atom_index(linker,"P")]
#    single_vector4 = LabF.generate_vector_of_interest(linker.get_O3PO5(), nextnuc.get_alpha(), [v_P, d2, d1])
#
#    # Add single_vector4 to , so that we define the location of O3'
#    v6 = LFT1.move_vector_to_loc(LFT1.return_normalized(single_vector4) * 1.6, v_P)
#    # O3' index
#    id_O3 = LFT2.retrieve_atom_index(nextnuc, "O3'")
#    # Get distance from nextnuc O3' to the position defined as O3'
#    O3_distance_from_origin = nextnuc.array[id_O3]
#    nextnuc_distance = v6 - O3_distance_from_origin
#    # Get the distance from the P_O3 now and add it to nextnuc_origin
#    #   which will move the entire molecule(nextnuc) to the position of O3'
#    nextnuc_loc = LFT1.move_vector_to_loc(nextnuc.array, nextnuc_distance)
#
##------------ Rotate nextnuc correctly along zeta
#    # O5' (v2), P (v3), O3' (v6), C3'
#    P_O3_C3 = 119.032 * (np.pi / 180)
#    zeta_dihr = nucleoside.get_zeta()
#
#    v2 = nucleotide[LFT2.retrieve_atom_index(nucleoside, "O5'") + linker.get_shape()]
#    single_vector5 = LabF.generate_vector_of_interest(P_O3_C3, zeta_dihr, [v6, v_P, d2])
#
#    # We have our vector, which is O3' -> C3', so now we rotate the the nextnuc onto it
#    # Get nextnuc's O3' atom to be in origin
#    distance_to_origin = nextnuc_loc[id_O3]
#    nextnuc_loc_tmp = LFT1.move_vector_to_origin(nextnuc_loc, distance_to_origin)
#
#    # rotate the vector O3' -> C3' onto single_vector5
#    id_C3 = LFT2.retrieve_atom_index(nextnuc, "C3'")
#    v7 = nextnuc_loc[id_C3]
#    O3_C3 = LFT1.return_normalized(v7 - v6)          # C3' - O3' gets the direction of O3' -> C3'
#    quaternion_zeta = LFT1.get_quaternion(single_vector5, O3_C3)
#
#    # Rotate nextnuc and move it into place
#    nextnuc_loc = LFT1.rotate_with_quaternion(quaternion_zeta, nextnuc_loc_tmp)
#    nextnuc_loc = nextnuc_loc + distance_to_origin
#
##------------ Rotate nextnuc correctly along epsilon
#    # P(v3) - O3'(v6) - C3'(v7) - C2' dihedral
#
#    # override the current C3' vector, since we changed its location
#    v7 = nextnuc_loc[id_C3]
#
#    # Get epsilon dihedral
#    epsilon_dihr = nucleoside.get_epsilon()
#    # Get angle O3' - C3' - C4'
#    O3_C3_C4 = 111.919 * (np.pi/180)
#
#    single_vector6 = LabF.generate_vector_of_interest(O3_C3_C4, epsilon_dihr, [v7, v6, v_P])
#
#    ## now that we have the vector, rotate the nucleoside appropriately.
#    # get C3' and use it as the distance of the atom of interest to the origin
#    C3_vector = nextnuc_loc[id_C3]
#    # get C4'
#    C4_vector = nextnuc_loc[LFT2.retrieve_atom_index(nextnuc, "C4'")]
#    nuc_at_OG = LFT1.move_vector_to_origin(nextnuc_loc, C3_vector)
#    nuc_vector_C3C4 = LFT1.return_normalized(C4_vector - C3_vector)                 # C4' - C3' gives C3' -> C4'
#
#    O3_C3 = LFT1.return_normalized(nuc_at_OG[id_C3] - nuc_at_OG[id_O3])             # needs to be in O3' -> C3' for some reason
#
#    # We need to rotate around the direction of the vector O3' -> C3'
#    quaternion_epsilon = LFT1.get_quaternion_custom_axis(single_vector6, nuc_vector_C3C4,rotation_axis = O3_C3)
#    rotated_nucatOG = LFT1.rotate_with_quaternion(quaternion_epsilon, nuc_at_OG)
#    nextnuc_loc = rotated_nucatOG + C3_vector


######################### CREATE THE PDB THAT GOES WITH ARRAY INPUTTED ####################
    leading_strand = (nextnuc_loc, leading_strand)
    LabF.create_PDB_from_matrix(np.vstack(leading_strand), nucleic_acid_list)


########################                    FINAL                    ###################### 
