import numpy as np


import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2
import labyrinth_func_tools3 as LFT3







""" This file is used to refactor some of the code from labyrinth_func.py, to make a more legible. """
def retrieve_index_of_best_conformation(stored_bb_distances, stored_bb_bools):

    # If this case happens, then check which one is closest to the defaulted size between the linker and the previous nucleotide
    if np.all(checkBool == True for checkBool in stored_bb_bools):
        diff_distances = np.zeros(shape=len(stored_bb_bools), dtype=float)
        for i in range(len(diff_distances)):
            diff_distances[i] = abs(stored_bb_distances[i] - 1.6)

        smallest_diff = diff_distances.min()
        index_smallest_dif = np.where(diff_distances == smallest_diff)[0]
        return int(index_smallest_dif)

    # Check if some are True. If there is only one True, just take that conformation. If there are multiple Trues, test them against eachother
    if np.any(checkBool == True for checkBool in stored_bb_bools):
        if np.count_nonzero(stored_bb_bools == True) == 1:
            only_true_conf = np.where(True == stored_bb_distances)[0]
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
