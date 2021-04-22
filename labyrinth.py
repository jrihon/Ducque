import numpy as np
import json
import sys
import os
import labyrinth_func
from matplotlib import pyplot as plt


""" Create dictionary of the filename and their molecule code """

codex_acidum_nucleicum = {
'DT': ['json/dna_thymidine.json', 'json/phosphate_linker.json']


}

def printexit(param):
    print(param)
    exit()


def Architecture(nucleic_acid):

    ## PARSE DATA FROM JSON
    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
    nucleo = labyrinth_func.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

    # Parse dictionary for the correct linker segment
    link = labyrinth_func.Desmos(codex_acidum_nucleicum[nucleic_acid][1])


    # Calculate the normal of C4'-C5'-O5' 
    v0 = nucleo.array[0] # O5'
    v1 = nucleo.array[1] # C5'

    C5O5 = (v1 - v0) * -1.0        # I-hat , where origin is equal to v1

    # Get vector of normalized magnitude 
    C5O5 = C5O5 / np.linalg.norm(C5O5)

    # angle in RADIANS
    _angleCOP = link.get_COP() * (np.pi / 180)
    # Generate the cone vector
    cone_vector = labyrinth_func.generate_cone_vector(_angleCOP)

    # Check if cone vector has the correct angles
    labyrinth_func.check_phi_angle_of_vector(cone_vector)

    # Rotate towards a direction, with a certain angle.
    rot_matrix = labyrinth_func.get_rM(C5O5)

    # Get the transformed cone of vectors
    rotated_vector = labyrinth_func.rotate_cone_vector(rot_matrix, cone_vector)

    # Check if all the vectors have been rotated correctly (OPTIONAL)
    labyrinth_func.check_phi_angle_of_vector(rotated_vector, axis=C5O5)

    # Calculate dihedrals with the cone
    range_of_dihedrals = labyrinth_func.praxeolitic_dihedralRANGE(nucleo.array, rotated_vector)

    # DIHEDRAL TO FIT THE VECTOR ON
    # interpolation
    y_interpolate = labyrinth_func.get_interpolated_dihedral(range_of_dihedrals, nucleo.get_beta())
    # generation
    single_vector = labyrinth_func.generate_rotate_single_vector(rot_matrix, y_interpolate, _angleCOP)

    # Check it's dihedral is correct (OPTIONAL)
    checkdihr = labyrinth_func.praxeolitic_dihedralSINGLE(nucleo.array, single_vector)
    print(nucleo.get_beta(), checkdihr)
    # Parse dictionary for the correct linker segment
    #link = labyrinth_func.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#
#    # vector of interest
#    ax.scatter(C5O5[0], C5O5[1], C5O5[2], color='gray')
#    #rotated cone
#    ax.scatter(rotated_vector[:,0], rotated_vector[:,1], rotated_vector[:,2], color='green', alpha=0.5)
#
#    #original cone
#    ax.scatter(cone_vector[:,0], cone_vector[:,1], cone_vector[:,2], color='r', alpha=0.5)
#
#    #origin
#    ax.scatter(0,0,0, color='black')
#
#    #rotated_vector
#    ax.scatter(single_vector[0], single_vector[1], single_vector[2], color='blue')
#
#    ax.set_zlim(-2,2)
#    ax.set_xlim(-2,2)
#    ax.set_ylim(-2,2)
#    ax.set_xlabel('x_axis')
#    ax.set_ylabel('y_axis')
#    ax.set_zlabel('z_axis')
#
#    plt.show()
