import numpy as np
import json
import sys
import os
import labyrinth_func as LabF
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
    nucleo = LabF.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

    # Parse dictionary for the correct linker segment
    link = LabF.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

    ## Make the required vectors and retrieve the data
    # Calculate the normal of C4'-C5'-O5' 
    v0 = nucleo.array[0] # O5'
    v1 = nucleo.array[1] # C5'
    C5O5 = (v1 - v0) * -1.0        # I-hat , where origin is equal to v1

    # Get vector of normalized magnitude 
    C5O5 = C5O5 / np.linalg.norm(C5O5)

    # angle of C5' - O5 ' - P in RADIANS
    _angleCOP = link.get_COP()

    # Generate the cone vector
    cone_vector1 = LabF.generate_cone_vector(_angleCOP)

    # Check if cone vector has the correct angles
    #LabF.check_phi_angle_of_vector(cone_vector)

    ## ROTATE
    # Rotate towards a direction, with a certain angle. Invert the vector, to get the correct angle for the cones
    # The angle will then be between vector( O5' - C5') and vector( O5' - P), with the same starting point.
    _quat1 = LabF.get_quaternion(C5O5 * -1.0)

    # Get the transformed cone of vectors
    rotated_vector1 = LabF.rotate_with_quaternion(_quat1, cone_vector1)

    # Check if all the vectors have been rotated correctly (OPTIONAL)
    LabF.check_phi_angle_of_vector(rotated_vector1, axis=C5O5 * -1.0)

    # Calculate dihedrals with the cone
    range_of_dihedrals1 = LabF.praxeolitic_dihedralRANGE([nucleo.array[0], nucleo.array[1], nucleo.array[4]], rotated_vector1)

    ## EXTRAPOLATE TO THE DIHEDRAL AND GENERATE THE NECESSARY VECTOR
    # interpolation
    theta_interpolate1 = LabF.get_interpolated_dihedral(range_of_dihedrals1, nucleo.get_beta())
    # generate
    single_vector1 = LabF.generate_and_rotate_single_vector_QUAT(theta_interpolate1, _angleCOP, _quat1)

    #(OPTIONAL) Check it's dihedral is correct
    #checkdihr = LabF.praxeolitic_dihedralSINGLE([nucleo.array[0], nucleo.array[1], nucleo.array[4]], single_vector1)

    # POSITION THE LINKER TO THE LOCATION THE NEW VECTOR POINTS AT
    # Move the linker to the correct location
    OP2_loc = LabF.position_linker(v0, single_vector1, link)

    #check_vector = (nucleo.array[0] - OP2_loc[0]) * -1.0
    #check_position = LabF.praxeolitic_dihedralSINGLE([ nucleo.array[0], nucleo.array[1], nucleo.array[4] ], check_vector)
#---------------------------------------------- SO FAR THE LINKER HAS BEEN POSITIONED, BUT IT NEEDS TO ROTATE  
    # The distance from P to the origin of the xyz system
    p_to_origin = OP2_loc[0]

    # move the linker to the origin, by positioning the posphorus at [0,0,0]
    OP2_origin = OP2_loc - p_to_origin
    p0 = OP2_origin[0]                 # P atom, with start at origin
    p1 = v0 - p_to_origin           # O5', with start at origin

    # Create vector on which we will rotate on vector(O5' - P) (when doing * -1.0 ; if not then P - O5') 
    O5P = (p1 - p0)
    # Normalise O5P
    O5P /= np.linalg.norm(O5P)

    # get angle of O5' - P - OP2 in RADIANS from the linker.json
    _angleOPO = link.get_OPO2()

    # Generate the cone vector, which gets generated around Z_axis
    cone_vector2 = LabF.generate_cone_vector(_angleOPO)
    LabF.check_phi_angle_of_vector(cone_vector2)

    # Rotate the cone towards a certain angle
    _quat2 = LabF.get_quaternion(O5P)

    # Get the transformed cone of vectors
    rotated_cone2 = LabF.rotate_with_quaternion(_quat2, cone_vector2)

    # Check if all vectors have been rotated correctly
    LabF.check_phi_angle_of_vector(rotated_cone2, axis=O5P)

    # Calculate the dihedrals with the cone
    range_of_dihedrals2 = LabF.praxeolitic_dihedralRANGE([OP2_loc[0], nucleo.array[0], nucleo.array[1]], rotated_cone2)
    ## DIHEDRAL TO FIT THE VECTOR ON
    # interpolation
    theta_interpolate2 = LabF.get_interpolated_dihedral(range_of_dihedrals2, link.get_OP2_dihedral())

    # generate and position the vector correctly
    single_vector2 = LabF.generate_and_rotate_single_vector_QUAT(theta_interpolate2, _angleOPO, _quat2)

    # I need P, O5', C5' to input. OP2 is the rotated vector
    checkdihr = LabF.praxeolitic_dihedralSINGLE([OP2_loc[0], nucleo.array[0], nucleo.array[1]], single_vector2)
    #printexit(checkdihr, link.get_OP2_dihedral())
        #-------------------------------------------------- rotate the phosphate group
    ## Rotate the phosphate to the single_vector
    # Define the vector that goes from P to OP2 and normalize it
    Phosp = OP2_loc[0]
    OP2 = OP2_loc[2]
    P_OP2 = (OP2 - Phosp)
    P_OP2 /= np.linalg.norm(P_OP2)

    # Get rotMatrix and align my phosphate group - vector with the single vector; matrix rotation
    _quat_phosphate = LabF.get_quaternion(single_vector2 , vector_to_rotate_from=P_OP2)
    #Rotate it
    OP2_loc2 = LabF.rotate_with_quaternion(_quat_phosphate, OP2_origin)

    # Move the rotated linker to the place of origin
    OP2_loc2 = OP2_loc2 + p_to_origin
    # DONE
    #------------------------------ROTATE THE LINKER AGAIN, BUT ON THE AXIS OF P_OP2, to get P_OP1

    ## Generate a cone vector that rotates on the axis of O5P again, but this time with the angle of P_OP1
    # Get the linker to the origin and extract the values of interest
    OP2_loc3 = OP2_loc2 - p_to_origin

    # the P_O2 vector ( OP2 - P resulteert in P -> OP2 )
    P_O2 = OP2_loc3[2] - OP2_loc3[0]
    P_O2 /= np.linalg.norm(P_O2)
    # the P_O1 vector
    P_O1 = OP2_loc3[1] - OP2_loc3[0]
    P_O1 /= np.linalg.norm(P_O1)
    # Calculate angle between these two and you will get the angle of OP1_P_OP2. 
    # The order here does not matter
    # the angle is automatically in radians
    angleLINKER = link.get_OPO1()
    #angle between the two oxygens
    #angleLINKER = np.arccos(np.dot(P_O1, P_O2)) * 180/np.pi

    # generate cone vector
    cone_vector3 = LabF.generate_cone_vector(angleLINKER)

    # generate the quaternion
    _quat3 = LabF.get_quaternion(O5P)
    #_quat3 = LabF.get_quaternion(P_O2)

    # Rotate the cone vector on top of P_O2
    rotated_cone3 = LabF.rotate_with_quaternion(_quat3, cone_vector3)
    #LabF.check_phi_angle_of_vector(rotated_cone3, axis=O5P)

    ## Interpolate the correct theta value of the dihedral we want
    # Calculate the dihedrals with the cone, here is P - O5' - C5' and then all the possible OP1 oxygens
    range_of_dihedrals3 = LabF.praxeolitic_dihedralRANGE([OP2_loc[0], nucleo.array[0], nucleo.array[1]], rotated_cone3)

    # interpolation
    theta_interpolate3 = LabF.get_interpolated_dihedral(range_of_dihedrals3, link.get_OP1_dihedral())

    # generate and position the vector correctly
    single_vector3 = LabF.generate_and_rotate_single_vector_QUAT(theta_interpolate3, angleLINKER, _quat3)

    # I need P, O5', C5' to input. OP2 is the rotated vector
    #checkdihr = LabF.praxeolitic_dihedralSINGLE([OP2_loc[0], nucleo.array[0], nucleo.array[1]], single_vector3)
    #checkangle = LabF.get_angle_for_rM(single_vector3, O5P) * (180/np.pi)
    #printexit(checkangle, link.get_OPO1() * (180/np.pi))
    #printexit(checkdihr, link.get_OP1_dihedral())
    #------------ Now that we have the vector of interest, we rotate OP_1 onto it

    # generate quaternion
    #_quat4 = LabF.get_quaternion(single_vector3, P_O1)
    _quat4 = LabF.get_quaternion_custom_axis(single_vector3, P_O1, P_O2 * -1.0)

    # rotate the vector
    OP_final = LabF.rotate_with_quaternion(_quat4, OP2_loc3)
    # Bring OP_final to the location of the phosphorus
    OP_final = OP_final + p_to_origin


######################### LINKER POSITIONED ###############################################


######################### POSITION THE NEXT NUCLEOTIDE ####################################




    LabF.create_PDB_from_matrix(np.vstack((OP_final, nucleo.array)), nucleo, link)
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
#    ax.scatter(single_vector[0], single_vector[1], single_vector[2], color='green')
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
