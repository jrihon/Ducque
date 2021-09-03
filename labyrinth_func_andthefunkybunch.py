import numpy as np

import labyrinth_func as LabF
import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2
import labyrinth_func_tools3 as LFT3



""" This file is used to refactor some of the code from labyrinth_func.py, to make a more legible. """




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

        multiple_true_conf = np.where(True == stored_distances)[0]
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



def assert_the_dihedral_of_interest(compl_nuc, compl_nuc_arr : np.ndarray, compl_linker, prev_nuc, complementary_strand : np.ndarray, index_compl, prev_linker) -> bool:
    # Atom Parsing List for the knowing which atoms to parse from the respective arrays ; for bond length and dihedral evaluation
    APL = LFT3.Atom_Parsing_List(compl_nuc, compl_linker, prev_nuc)
    dihedral = compl_nuc.get_dihedral("alpha")
    angle = compl_nuc.get_angle("alpha")

    # Parse the last atom needed from the complementary strand. NB : index_compl has been set to the correct value in the function ' assert_possible ... _and_fit() '.
    id_compl_strand = LFT2.retrieve_atom_index(prev_nuc, APL[3]) + index_compl + prev_linker.mol_length
    v_compl_strand = complementary_strand[id_compl_strand]

    # Parse the other atoms needed from the current nucleotide that is being fitted
    id_v0 = LFT2.retrieve_atom_index(compl_nuc, APL[0]) + compl_linker.mol_length
    id_v1 = LFT2.retrieve_atom_index(compl_nuc, APL[1]) + compl_linker.mol_length
    id_v2 = LFT2.retrieve_atom_index(compl_linker, APL[2])

    v0 = compl_nuc_arr[id_v0]
    v1 = compl_nuc_arr[id_v1]
    v2 = compl_nuc_arr[id_v2]

    # Calculate the dihedral. Let's say the dihedral should not deviate more than 25 degrees?
    calculated_dihr = LFT1.dihedral_single(v0, v1, v2, v_compl_strand)

    # Create booleans to assert the calculated values and see if they are within ranges of concordance
    bool_dihedral = LFT1.assert_dihedral_of_nucleotide(calculated_dihr, dihedral)

    return bool_dihedral


def capping_retrieve_atomnames(list_of_leading_sequence : list, list_of_complementary_sequence : list, backbone_dict : dict, nuc_dict : dict) -> list:
    """ This should return the names of the capping atoms.
            i.e. ; if the H binds to O5', the atom name of the H should be HO5' 

        Parse the backbone atoms from the respective dictionaries over at labyrinth_func_tools3.py """
    # Initiate a list with filler values
    list_lead_atom = [0, 0]
    list_compl_atom = [0, 0]
    index_list = [-1, 0]
    for atom in range(len(index_list)):
        # Leading strand
        lead_nucleoside_fname = nuc_dict[list_of_leading_sequence[index_list[atom]]][0]
        atom_lead = LabF.Nucleoside(lead_nucleoside_fname)
        # Complementary strand
        compl_nucleoside_fname = nuc_dict[list_of_complementary_sequence[index_list[atom]]][0]
        atom_compl = LabF.Nucleoside(compl_nucleoside_fname)

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
        lead_atom = LabF.Nucleoside(lead_nucleoside_fname)
        # Complementary strand
        compl_nucleoside_fname = nuc_dict[list_of_complementary_sequence[index_list[i]]][0]
        compl_atom = LabF.Nucleoside(compl_nucleoside_fname)

        backbone_lead = backbone_dict[lead_atom.get_nucleic_acid_code()]
        backbone_compl = backbone_dict[compl_atom.get_nucleic_acid_code()]

        # Now retrieve the index of the arrays for which we want to parse the atoms
        if index_list[i] == 0 :
            # Parse the atom of interest from the backbone
            for j in [-1, -2, -3]:
                lead_atom_list[i][j] = backbone_lead[j]
                compl_atom_list[i][j] = backbone_compl[j]

            lead_atom_indexes = LFT2.retrieve_atom_index_MULTIPLE(lead_atom, lead_atom_list[i])
            compl_atom_indexes = LFT2.retrieve_atom_index_MULTIPLE(compl_atom, compl_atom_list[i])

            # Now retrieve the coordinates from the array
            l0 = leading_array[lead_atom_indexes[-3]]
            l1 = leading_array[lead_atom_indexes[-2]]
            l2 = leading_array[lead_atom_indexes[-1]]

            c0 = complementary_array[compl_atom_indexes[-3]]
            c1 = complementary_array[compl_atom_indexes[-2]]
            c2 = complementary_array[compl_atom_indexes[-1]]

            # Now calculate the would be dihedral and retrieve the H position in xyz
            vector_lead = LabF.generate_vector_of_interest(angle, dihedral, [l2, l1, l0])
            vector_compl = LabF.generate_vector_of_interest(angle, dihedral, [c2, c1, c0])

            # Correct the length, add the vector to the array and store it in H_vectors
            H_vectors[0] = (LFT1.return_normalized(vector_lead) * 1.2) + l2
            H_vectors[2] = (LFT1.return_normalized(vector_compl) * 1.2) + c2

        elif index_list[i] == -1 :
            for j in [0, 1, 2]:
                lead_atom_list[i][j] = backbone_lead[j]
                compl_atom_list[i][j] = backbone_compl[j]

            # Create the index counter for both lead and compl strand
            lead_index = leading_array.shape[0] - lead_atom.mol_length
            compl_index = complementary_array.shape[0] - compl_atom.mol_length

            lead_atom_indexes = LFT2.retrieve_atom_index_MULTIPLE(lead_atom, lead_atom_list[i], lead_index)
            compl_atom_indexes = LFT2.retrieve_atom_index_MULTIPLE(compl_atom, compl_atom_list[i], compl_index)

            # Now retrieve the coordinates from the array
            l0 = leading_array[lead_atom_indexes[0]]
            l1 = leading_array[lead_atom_indexes[1]]
            l2 = leading_array[lead_atom_indexes[2]]

            c0 = complementary_array[compl_atom_indexes[0]]
            c1 = complementary_array[compl_atom_indexes[1]]
            c2 = complementary_array[compl_atom_indexes[2]]

            # Now calculate the would be dihedral and retrieve the H position in xyz
            vector_lead = LabF.generate_vector_of_interest(angle, dihedral, [l0, l1, l2])
            vector_compl = LabF.generate_vector_of_interest(angle, dihedral, [c0, c1, c2])

            # Correct the length, add the vector to the array and store it in H_vectors
            H_vectors[1] = (LFT1.return_normalized(vector_lead) * 1.2) + l0
            H_vectors[3] = (LFT1.return_normalized(vector_compl) * 1.2) + c0


    return H_vectors

