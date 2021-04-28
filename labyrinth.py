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

def printexit(*arg):
    """Function to print out any numbers of parameters and then exit the software """
    print(arg)

    sys.exit(0)


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

    # angle of C5' - O5 ' - P in RADIANS
    _angleCOP = link.get_COP()
    # Generate the cone vector
    cone_vector = labyrinth_func.generate_cone_vector(_angleCOP)

    # Check if cone vector has the correct angles
    #labyrinth_func.check_phi_angle_of_vector(cone_vector)

    # Rotate towards a direction, with a certain angle. Invert the vector, to get the correct angle for the cones
    # The angle will then be between vector( O5' - C5') and vector( O5' - P), with the same starting point.
    rot_matrix = labyrinth_func.get_rM(C5O5 * -1.0)

    # Get the transformed cone of vectors
    rotated_vector = labyrinth_func.rotate_cone_vector(rot_matrix, cone_vector)

    # Check if all the vectors have been rotated correctly (OPTIONAL)
    labyrinth_func.check_phi_angle_of_vector(rotated_vector, axis=C5O5 * -1.0)

    # Calculate dihedrals with the cone
    range_of_dihedrals = labyrinth_func.praxeolitic_dihedralRANGE([nucleo.array[0], nucleo.array[1], nucleo.array[4]], rotated_vector)
    # DIHEDRAL TO FIT THE VECTOR ON
    # interpolation
    theta_interpolate = labyrinth_func.get_interpolated_dihedral(range_of_dihedrals, nucleo.get_beta())
    # generate
    single_vector = labyrinth_func.generate_rotate_single_vector(theta_interpolate, _angleCOP, rot_matrix)

    # Check it's dihedral is correct (OPTIONAL)
    checkdihr = labyrinth_func.praxeolitic_dihedralSINGLE([nucleo.array[0], nucleo.array[1], nucleo.array[4]], single_vector)

    # Move the linker to the correct location
    OP2_loc = labyrinth_func.position_linker(v0, single_vector, link)

    check_vector = (nucleo.array[0] - OP2_loc[0]) * -1.0

    check_position = labyrinth_func.praxeolitic_dihedralSINGLE([ nucleo.array[0], nucleo.array[1], nucleo.array[4] ], check_vector)
#---------------------------------------------- SO FAR SO GOOD  
    # The first index of the linker is the phosphorus
    p0 = OP2_loc[0]     # P
    p1 = v0             # O5'
    # Create vector on which we will rotate on vector(O5' - P)
    O5P = (p1 - p0) * -1.0

    # Normalise O5P
    O5P /= np.linalg.norm(O5P)

    # get angle of O5' - P - OP2 in RADIANS
    _angleOPO = link.get_OPO()

    # Generate the cone vector, which gets generated around Z_axis
    cone_vector = labyrinth_func.generate_cone_vector(_angleOPO)
    labyrinth_func.check_phi_angle_of_vector(cone_vector)

    # Rotate the cone towards a certain angle
    rot_matrix = labyrinth_func.get_rM(O5P * -1.0)

    # Get the transformed cone of vectors
    rotated_cone = labyrinth_func.rotate_cone_vector(rot_matrix, cone_vector)

    # Check if all vectors have been rotated correctly
    labyrinth_func.check_phi_angle_of_vector(rotated_cone, axis=O5P * -1.0)

    # Calculate the dihedrals with the cone
    range_of_dihedrals = labyrinth_func.praxeolitic_dihedralRANGE([OP2_loc[0], nucleo.array[0], nucleo.array[1]], rotated_cone)
    #DIHEDRAL TO FIT THE VECTOR ON
    # interpolation
    # I need P, O5', C5' to input. OP2 is the rotated vector
    y_interpolate = labyrinth_func.get_interpolated_dihedral(range_of_dihedrals, link.get_OP2_dihedral())
    # generate

    single_vector = labyrinth_func.generate_rotate_single_vector(y_interpolate, _angleOPO, rot_matrix)
    # Check if it's done correctly
    checkdihr = labyrinth_func.praxeolitic_dihedralSINGLE([OP2_loc[0], nucleo.array[0], nucleo.array[1]], single_vector)
    #printexit(checkdihr, link.get_OP2_dihedral())
#-------------------------------------------------- rotate the phosphate group
    ## Rotate the phosphate to the single_vector
    # Get the phosphorus linker to the origin
    p_to_origin = OP2_loc[0]

    #print(np.degrees(np.arccos(np.dot(single_vector, (OP2_loc[0] - nucleo.array[0])*-1.0)))) ; exit()
    OP2_loc -= p_to_origin
    # Define the vector that goes from P to OP2
    Phosp = OP2_loc[0]
    OP2 = OP2_loc[2]
    P_OP2 = (OP2 - Phosp)
    P_OP2 /= np.linalg.norm(P_OP2)

#
#    # Get rotMatrix and align my phosphate group - vector with the single vector; matrix rotation
#    rM_phosphate = labyrinth_func.get_rM(single_vector, vector_to_rotate_from=P_OP2)
#
#    OP2_loc = labyrinth_func.rotate_single_vector(rM_phosphate, OP2_loc)
#
#
#    # DONE

    labyrinth_func.create_PDB_from_matrix(np.vstack((OP2_loc, nucleo.array)), nucleo, link) 
#    # Turn the segment to the correct dihedral
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#
#    # vector of interest
#    #ax.scatter(C5O5[0], C5O5[1], C5O5[2], color='gray')
#    #rotated cone
#    #ax.scatter(rotated_vector[:,0], rotated_vector[:,1], rotated_vector[:,2], color='green', alpha=0.5)
#
#    #original cone
#    #ax.scatter(cone_vector[:,0], cone_vector[:,1], cone_vector[:,2], color='r', alpha=0.5)
#
#    #origin
#    ax.scatter(0,0,0, color='black')
#
#    #rotated_vector
#    ax.scatter(single_vector[0], single_vector[1], single_vector[2], color='orange')
#    ax.scatter(nucleo.array[:,0], nucleo.array[:,1],nucleo.array[:,2], color='blue')
#    ax.scatter(OP2_loc[:,0], OP2_loc[:,1], OP2_loc[:,2], color='orange')
#
##    ax.set_zlim(-2,2)
##    ax.set_xlim(-2,2)
##    ax.set_ylim(-2,2)
#    ax.set_xlabel('x_axis')
#    ax.set_ylabel('y_axis')
#    ax.set_zlabel('z_axis')

#    plt.show()
