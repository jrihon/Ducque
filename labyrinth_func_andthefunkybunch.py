import numpy as np

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
