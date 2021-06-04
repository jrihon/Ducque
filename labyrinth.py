import sys, os, json
import numpy as np
import labyrinth_func as LabF


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
    # Get vector of normalized magnitude 
    C5O5 = LabF.return_normalized((v1 - v0) * -1.0)

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
    O5P = LabF.return_normalized(p1 - p0)

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
    P_OP2 = LabF.return_normalized(OP2 - Phosp)

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
    P_O2 = LabF.return_normalized(OP2_loc3[2] - OP2_loc3[0])
    # the P_O1 vector
    P_O1 = LabF.return_normalized(OP2_loc3[1] - OP2_loc3[0])
    # Get angle
    angleLINKER = link.get_OPO1()

    # generate cone vector
    cone_vector3 = LabF.generate_cone_vector(angleLINKER)

    # generate the quaternion
    _quat3 = LabF.get_quaternion(O5P)

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
    _quat4 = LabF.get_quaternion_custom_axis(single_vector3, P_O1, P_O2 * -1.0)

    # rotate the vector
    OP_final = LabF.rotate_with_quaternion(_quat4, OP2_loc3)
    # Bring OP_final to the location of the phosphorus
    OP_final = OP_final + p_to_origin


######################### LINKER POSITIONED ###############################################


######################### POSITION THE NEXT NUCLEOTIDE ####################################
# Start with a the stacked array, that makes it easier
    nucleotide = np.vstack((OP_final, nucleo.array))

    # import the next nucleoside
    nextnuc = LabF.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

#--------------- We start with creating the position for the O3'
    # We create the two vectors, of the three, that make up the dihedral. Then we generate a set of vectors and look for the third vector.
    C5O5_sec = LabF.return_normalized(nucleotide[3] - nucleotide[4])      # O5' - C5' makes C5' -> O5' direction
    O5P_sec = LabF.return_normalized(nucleotide[0] - nucleotide[3])       # P - O5' makes O5' -> P direction

    # We need the angle
    _angleO3PO5 = link.get_O3PO5()

    # We need the dihedral angle of interest
    dihr_alpha = nextnuc.get_alpha()

    # We need the dihedral angle, so we start off by generating a cone
    cone_vector_next1 = LabF.generate_cone_vector(_angleO3PO5)

    # Rotate the cone onto the vector of interest. We invert to get the angle between O5' - P - O3'
    _quat_next1 = LabF.get_quaternion(O5P_sec * -1.0)

    # Turn the cone_vector
    rotated_cone_next1 = LabF.rotate_with_quaternion(_quat_next1, cone_vector_next1)

    # Calculate the set of dihedral angles. Input the vectors in the reverse direction.
    # if ( C5' - O5' - P - O3'), then ([P, O5', C5'], conevector O3')
    range_of_dihedrals_next1 = LabF.praxeolitic_dihedralRANGE([nucleotide[0], nucleotide[3], nucleotide[4]], rotated_cone_next1)

    # Get the interpolated dihedral angle
    theta_interpolate_next1 = LabF.get_interpolated_dihedral(range_of_dihedrals_next1, dihr_alpha)

    # generate and position the vector correctly
    single_vector3 = LabF.generate_and_rotate_single_vector_QUAT(theta_interpolate_next1, _angleO3PO5, _quat_next1)

    #------------------- Ok... now we have the vector
    #   We create the vector P -> 03' and turn this vector onto single_vector3
    #   Then, we we turn the entire nextnuc on top of this vector

    # add single_vector3 to the end of phosphate linker, so that we position the location of O3'
    P_O3 = LabF.resultant_vector_addition(LabF.return_normalized(single_vector3) * 1.6, nucleotide[0])
    # Get distance from nextnuc O3' to the position defined as O3'
    nextnuc_distance = P_O3 - nextnuc.array[28]
    # Get the distance from the P_O3 now and add it to nextnuc_origin
    #   which will move the entire molecule(nextnuc) to the position of O3'
    nextnuc_loc = LabF.resultant_vector_addition(nextnuc.array, nextnuc_distance)

    #### DONE

#### OK, NOW WE GO FOR THE ZETA DIHEDRAL
        # O5', P, O3', C3' are what the zeta dihedral consists of
    O5_P3 = LabF.return_normalized(nucleotide[0] - nucleotide[3])       # P - O5' gets the direction of O5' -> P                    
    P_O33 = LabF.return_normalized(P_O3 - nucleotide[0])                # O3' - P gets the direction of P -> O3'
        # Get zeta dihedral value
    zeta_dihr = nucleo.get_zeta()
        # Get angle P_O3'_C3'
    P_O3_C3 = 119.032 * (np.pi / 180)
        # generate cone
    cone_vector_next2 = LabF.generate_cone_vector(P_O3_C3)
        # Get quaternion
    _quat_next2 = LabF.get_quaternion(P_O33 * -1.0)
        # Turn cone
    rotated_cone_next2 = LabF.rotate_with_quaternion(_quat_next2, cone_vector_next2)
        # interpolate dihedral. Zeta is O5' - P - O3' - C3', then [O3', P, O5']
    range_of_dihedrals_next2 = LabF.praxeolitic_dihedralRANGE([P_O3, nucleotide[0],nucleotide[3]] , rotated_cone_next2)
        # Get the interpolated dihedral angle
    theta_interpolate_next2 = LabF.get_interpolated_dihedral(range_of_dihedrals_next2, zeta_dihr)
        #LabF.check_phi_angle_of_vector(rotated_cone_next2, P_O33 * -1.0)
    single_vector4 = LabF.generate_and_rotate_single_vector_QUAT(theta_interpolate_next2, P_O3_C3, _quat_next2)

    # We have our vector, which is O3' -> C3', so now we rotate the the nextnuc onto it
    # Get nextnuc's O3' atom to be in origin
    dist_to_OG = nextnuc_loc[28]
    nextnuc_locAtOG = LabF.move_vector_to_origin(nextnuc_loc, dist_to_OG)

    # Rotate nextnuc_locAtOg onto the single_vector4
    O3_C3 = LabF.return_normalized(nextnuc_locAtOG[23] - nextnuc_locAtOG[28])          # C3' - O3' gets the direction of O3' -> C3'
    _quat_zeta = LabF.get_quaternion(single_vector4, O3_C3)

    # Rotate nextnuc and move it into place
    nextnuc_loc = LabF.rotate_with_quaternion(_quat_zeta, nextnuc_locAtOG)
    nextnuc_loc = nextnuc_loc + dist_to_OG

### OK, NOW WE GO FOR THE EPSILON DIHEDRAL
    P_O34 = LabF.return_normalized(nextnuc_loc[28] - nucleotide[0])            # O3' - P returns a vector P -> O3' 
    O3_C34 = LabF.return_normalized(nextnuc_loc[23] - nextnuc_loc[28])           # C3' - O3' returns a vector O3' -> C3'
        #get epsilon dihedral
    epsilon_dihr = nucleo.get_epsilon()
    print(epsilon_dihr)
        #get angle O3' - C3' - C4'
    O3_C3_C4 = 111.919 * (np.pi/180)
    increment_dihedral = 1 * (np.pi/180)
        #generate cone
    cone_vector_next3 = LabF.generate_cone_vector(O3_C3_C4)
        # get quaternion
    _quat_next3 = LabF.get_quaternion(O3_C34 * -1.0)
        # Turn cone
    rotated_cone_next3 = LabF.rotate_with_quaternion(_quat_next3, cone_vector_next3)
        # generate range of dihedral. Epsilon is P - O3' - C3' - C4', se reverse it [C3', O3', P]
    range_of_dihedrals_next3 = LabF.praxeolitic_dihedralRANGE([nextnuc_loc[23], nextnuc_loc[28], nucleotide[0]], rotated_cone_next3)
        # interpolate said range of dihedrals
    theta_interpolate_next3 = LabF.get_interpolated_dihedral(range_of_dihedrals_next3, epsilon_dihr)
        # generate single vector
    single_vector_next3 = LabF.generate_and_rotate_single_vector_QUAT(theta_interpolate_next3, O3_C3_C4, _quat_next3)
        #LabF.check_phi_angle_of_vector(rotated_cone_next3, O3_C34 * -1.0)

    ## now that we have the vector, rotate the nucleoside appropriately.
    # get C3' to the origin
    C3_from_OG = nextnuc_loc[23]
    nuc_at_OG = LabF.move_vector_to_origin(nextnuc_loc, C3_from_OG)
    nuc_vector_C3C4 = LabF.return_normalized(nuc_at_OG[4] - nuc_at_OG[23])              # C4' - C3' gives C3' -> C4'
    axis_of_rot = LabF.return_normalized(nuc_at_OG[28] - nuc_at_OG[23]) * -1.0                # needs to be in C3' -> O3'
    # We need to rotate around the direction of the vector C3' -> O3' (where we do O3_C34 * -1.0)
    _quat_next4 = LabF.get_quaternion_custom_axis(single_vector_next3, nuc_vector_C3C4, axis_of_rot)

    # turn the nucleoside and then put it back
    rotated_nucatOG = LabF.rotate_with_quaternion(_quat_next4, nuc_at_OG)
    nextnuc_loc = rotated_nucatOG + C3_from_OG




######################### CREATE THE PDB THAT GOES WITH ARRAY INPUTTED ####################
    stacked_array = (nextnuc_loc, OP_final, nucleo.array)
    LabF.create_PDB_from_matrix(np.vstack(stacked_array), nucleo, link)


########################                    FINAL                    ###################### 
