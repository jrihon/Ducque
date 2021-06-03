"""
The python file that organises the labyrinth_func.py better
"""
import numpy as np


def check_slope_of_array(arr : np.array) -> str:
    """ In labyrinth_func.py there is a function that retrieves the interpolated dihedral angle
    But it works on whether or not the list is ascending or descending
    That's what we need to figure out here now and return this """

    slope_list = []
    for i in range(len(arr)):
        if i == (len(arr) - 1):
            if arr[i] > arr[0]:
                slope_list.append("D")
            else:
                slope_list.append("A")

        elif arr[i] > arr[i+1]:
            slope_list.append("D")
        else:
            slope_list.append("A")

    descending = slope_list.count("D")
    ascending = slope_list.count("A")

    if ascending > descending:
        return "ASCENDING"
    else:
        return "DESCENDING"
