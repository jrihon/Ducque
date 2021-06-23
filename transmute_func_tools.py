import sys


""" Contains all the dictionaries and the tools to help the Transmutation work properly.    """

# Identify the nucleic acid chemistry
nucleoside_dict = {
        "DNA" : "Deoxyribonucleic acid",
        "RNA" : "Ribonucleic acid",
        }

linker_dict = {
        "DNA" : "Phosphate",
        "RNA" : "Phosphate"
        }


# Used for the filename and to fill out the identity dict in the json file
base_dict = {
        "A" : "Adenosine",
        "C" : "Cytidine",
        "G" : "Guanosine",
        "T" : "Thymidine",
        "U" : "Uracil",
        }


# The general backbone, whether it is a purine or a pyrimidine. This is a nested dictionary.
# Used to parse the indices of the vectors in the molecule's array
dihedral_dict = {
        "DNA" : {"backbone" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
                 "purine" : ["O4'", "C1'", "N9", "C4"],
                 "pyrimidine" : ["O4'", "C1'","N1", "C2"]
                 },
        "RNA" : {"backbone" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
                 "purine" : ["O4'", "C1'", "N9", "C4"],
                 "pyrimidine" : ["O4'", "C1'","N1", "C2"]
                 },
        }


# The general backbone, whether it is a purine or a pyrimidine. This is a nested dictionary.
# Used to parse the indices of the vectors in the molecule's array
angle_dict = {
        "DNA" : {"backbone" : [""],
                 "purine" : ["C1'", "N9", "C4"],
                 "pyrimidine" : ["C1'","N1", "C2"],
                 },
        "RNA" : {"backbone" : [""],
                 "purine" : ["C1'", "N9", "C4"],
                 "pyrimidine" : ["C1'","N1", "C2"],
                 }
        }


def get_base_type(base : str) -> str:
    """ Retrieve the type of base we will calculate with"""

    if base == "A":
        return "purine"

    if base == "T":
        return "purine"

    if base == "C":
        return "pyrimidine"

    if base == "G":
        return "pyrimidine"

    if base == "U":
        return "pyrimidine"

    # If you've reached this part, then something has gone wrong and it is not clear what the base is
    if True:
        print("It is not clear what the base is. Please check the Residue Name column in the prompted pdb file.\n"
                + "NB: The last character of the string should end as the identifier of one of the five canonical bases.")
        sys.exit(0)

#------------------------------------- NEVER USED THESE FUNCTIONS -------------------------------------#
#import numpy as np
#def get_indices_of_atoms(atom_names : list, atom_list : np.array) -> list:
#    """ retrieve the index of each respective atom in the atom list """
#    list_of_indices = []
#    for atom in atom_names:
#        ind = atom_list.index(atom)
#        list_of_indices.append(ind)
#
#    return list_of_indices
#
#def calculate_dihedral(atom_array : np.array, indices : list) -> float:
#
#    # Suppose the dihedral angle goes as follows : C4' - C5' - O5' - P
#
#    v2 = atom_array[indices[2]]     # O5'
#    v1 = atom_array[indices[1]]     # C5'
#    v0 = atom_array[indices[0]]     # C4'
#    v3 = atom_array[indices[4]]     # Supposed P
#
#    b0 = (v1 - v0) * -1.0           # C4' - C5' 
#    b1 = v2 - v1                    # C5' - O5'
#    b2 = v3 - v2                    # O5' - P
#
#    # normalize b1 so that it does not influence magnitude of vector rejections that come next
#    b1 /= LA.norm(b1)
#
#        # vector rejections
#    # v = projection of b0 onto plane perpendicular to b1
#    #   = b0 minus component that aligns with b1
#    # w = projection of b2 onto plane perpendicular to b1
#    #   = b2 minus component that aligns with b1
#    v = b0 - np.dot(b0, b1)*b1
#    w = b2 - np.dot(b2, b1)*b1
#
#    # angle between v and w in a plane is the torsion angle
#    # v and w may not be normalized but that's fine since tan is y/x
#    x = np.dot(v, w)
#    y = np.dot(np.cross(b1, v), w)
#
#    return np.degrees(np.arctan2(y, x))
#
#
#def calculate_angle(atom_array : np.array, indices : list) -> float:
#    """ The scalar product (dot product) to get the cosine angle.
#        Here we do the arccos, the get the angle immediately. """
#    a0 = atom_array[indices[0]]
#    a1 = atom_array[indices[1]]
#    a2 = atom_array[indices[2]]
#
#    c1 = a1 - a0
#    c2 = (a2 - a1) * -1.0
#
#    return np.arccos(np.dot(c1, c2))
#
