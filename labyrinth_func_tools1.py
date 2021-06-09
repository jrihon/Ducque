import numpy as np
import pandas as pd
import numpy.linalg as LA
from scipy.spatial.transform import Rotation as R

from typing import Tuple
import json
import labyrinth_func_tools2 as LFT2

"""
labyrinth_func_tools1 contains all the mathematical functions required, in python, to calculate for the rotations.
"""
# FUNCTIONS

def return_normalized(vector : np.ndarray) -> np.ndarray:
    """ returns a normalized vector """
    return vector / LA.norm(vector)


def generate_cone_vector(phi_angle : float) -> np.ndarray:
    """ Phi is the angle the vector makes with the Z-axis

    Graphing Spherical Coordinates in GeoGebra 3D (Part 2): A Cone about z-axis
    https://www.youtube.com/watch?v=Wl3Z3AqfI6c

    Generating uniform unit random vectors in R^n
    phi needs to be in RADIANS """

    rho = np.full(18, 1)
    theta = np.linspace(0, 2*np.pi, num=18, endpoint=False)
    phi = np.full(18, phi_angle)

    return np.array([ rho * np.cos(theta) * np.sin(phi),
                      rho * np.sin(theta) * np.sin(phi),
                      rho * np.cos(phi)]
                    ).T


def generate_and_rotate_single_vector_QUAT(interpolated_theta_angle : float, phi : float, quaternion) -> np.ndarray:
    """ Generate a single vector with the correct angle.
        Then rotate said angle to the correct orientation using the previously used quaternion.
        As always, phi needs to be in RADIANS"""

    single_vector = np.array([ 1.0 * np.cos(interpolated_theta_angle) * np.sin(phi),
                               1.0 * np.sin(interpolated_theta_angle) * np.sin(phi),
                               1.0 * np.cos(phi)]
                            ).T

    rotated_vector = quaternion.apply(single_vector)

    return rotated_vector


def get_direction_for_rM(from_vector : np.ndarray, vector_to_rotate_onto : np.ndarray) -> np.ndarray:
    """ cross product with the cone's axis to get the direction of the rotation axis
        Get the direction where vectorA rotates onto vectorB ; u = vectora X vectorb
        Let's normalize the direction """

    return np.cross(from_vector, vector_to_rotate_onto) / LA.norm(np.cross(from_vector, vector_to_rotate_onto))


def get_angle_for_rM(from_vector : np.ndarray, vector_to_rotate_onto : np.ndarray) -> float:
    """The scalar product (dot product) to get the cosine angle.
        Here we do the arccos, the get the angle immediately. """
    return np.arccos(np.dot(from_vector, vector_to_rotate_onto))


def get_quaternion(vector_to_rotate_onto : np.ndarray, vector_to_rotate_from : np.array = np.array([0,0,1] )):
    """
    Quaternion mathematics using the scipy.spatial.transform.Rotation library

    Create the quaternion that is associated with the angle and axis of rotation

    Apply it to the vector you want to rotate

    """
    # axis has been normalised
    axis = get_direction_for_rM(vector_to_rotate_from, vector_to_rotate_onto)      # DIRECTION

    # angle already in radians 
    theta = get_angle_for_rM(vector_to_rotate_from, vector_to_rotate_onto)         # ANGLE

    # Create the quaternion
    qx = axis[0] * np.sin(theta/2)
    qy = axis[1] * np.sin(theta/2)
    qz = axis[2] * np.sin(theta/2)
    qw = np.cos(theta/2)

    quaternion = R.from_quat([qx, qy, qz, qw])

    return quaternion


def get_quaternion_custom_axis(vector_to_rotate_onto : np.ndarray, vector_to_rotate_from : np.ndarray, rotation_axis : np.ndarray):
    """ Generate quaternion for when you already have the axis of rotation"""

    # normalise the axis
    axis = rotation_axis / LA.norm(rotation_axis)                            # DIRECTION
    # angle already in radians 
    theta = get_angle_for_rM(vector_to_rotate_from, vector_to_rotate_onto)       # ANGLE

    # Create the quaternion
    qx = axis[0] * np.sin(theta/2)
    qy = axis[1] * np.sin(theta/2)
    qz = axis[2] * np.sin(theta/2)
    qw = np.cos(theta/2)

    quaternion = R.from_quat([qx, qy, qz, qw])

    return quaternion


def rotate_with_quaternion(quaternion, vector : np.ndarray) -> np.ndarray:
    """ Vector rotation through quaternion mathematics"""

    return quaternion.apply(vector)


def check_phi_angle_of_vector(vectors : np.ndarray, axis : np.array = np.array([0,0,1])) -> None:
    """ The dotproduct determines the angle of the vector with a given axis/vector """

    ### to check of the angle is correct, you should introduce an if statement
    #   in the labyrinth code so that if the cone vector is not placed correctly
    #   that you just invert the vector you place it on (so vector * -1.0)
    ###


    # if a vector is given and it is not k^-hat
    if list(axis) != [0,0,1]:
        arr = np.empty(shape=(vectors.shape[0]))

        for i in range(vectors.shape[0]):
            arr[i] = np.arccos(np.dot(vectors[i], axis))

        #print(np.degrees(arr))
        arr = np.around(arr, decimals=10)                   # Round of values at the end to avoid errors on equality at the 17th decimal place
        if list(arr).count(arr[0]) != len(arr):
            print('Not all cone angles are correctly generated')
            sys.exit(0)

        return

    # If no vector is given, default to k^-hat; Z_axis)
    arr = np.arccos(np.dot(vectors, axis))
    #print(np.degrees(arr))
    arr = np.around(arr, decimals=10)
    if not list(arr).count(arr[0]) == len(arr):

        print('Not all cone angles are correctly generated')
        sys.exit(0)


def dihedral_array(atoms_in_dihr : np.ndarray, cone_vector : np.ndarray) -> np.ndarray :

    """
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    Praxeolitic dihedral calculation. as of now, the most efficient way to calculate for a dihedral angle.
    """

    v2 = atoms_in_dihr[0] # O5'
    v1 = atoms_in_dihr[1] # C5'
    v0 = atoms_in_dihr[2] # C4'
    v3 = cone_vector    # P , this contains all the possible P atoms

    b0 = -1.0 * (v1 - v0)   #C4-C5
    b1 = v2 - v1            #C5-O5
    b2 = v3                 #O5-P

    # normalize b1 so that it does not influence magnitude of vector rejections that come next
    b1 /= LA.norm(b1)

        # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    dihedrals = np.empty(b2.shape[0])

    for i in range(b2.shape[0]):
        v = b0 - np.dot(b0, b1)*b1
        w = b2[i] - np.dot(b2[i], b1)*b1

        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        dihedrals[i] = np.degrees(np.arctan2(y, x))

    return dihedrals


def dihedral_single(atoms_in_dihrjson_array : np.ndarray, single_vector : np.ndarray) -> np.ndarray:

    v2 = atoms_in_dihr[0] # O5'
    v1 = atoms_in_dihr[1] # C5'
    v0 = atoms_in_dihr[2] # C4'
    v3 = single_vector # Supposed P

    b0 = -1.0 * (v1 - v0)   #C4-C5 
    b1 = v2 - v1            #C5-O5
    b2 = v3                 #O5-P

    # normalize b1 so that it does not influence magnitude of vector rejections that come next
    b1 /= LA.norm(b1)

        # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return np.degrees(np.arctan2(y, x))


def interpolate_dihedrals(tuple_dihr : Tuple[Tuple[float, float], Tuple[float, float]], angle_dihr : float) -> float:
    """
    y = y1 + [ (x - x1) / (x2 - x1) * (y2 - y1) ]

    tuple_dihr = ( (rad1, degrees1), (rad2, degrees2) )

    The y's are the theta angles, indexed by the keys of the dictionary being parsed
    """
    # Print out the tuple dihedral and the angle. Just to see how it comes out
    #print(tuple_dihr, angle_dihr)

    # Calculate the interpolated angle
    y2 = tuple_dihr[0][0]           # radians in the cone generator function after dihedral of interest
    y1 = tuple_dihr[1][0]           # radians in the cone generator functionbefore dihedral of interest
    x2 = tuple_dihr[0][1]           # value of degrees in the cone generator function after ..
    x1 = tuple_dihr[1][1]           # value of degrees in the cone generator function before ...
    x = angle_dihr                  # dihedral angle of interest, used to extrapolate to the exact radians value in the cone generator

    return y1 + ( (x - x1) / (x2 - x1) * (y2 - y1) )


def get_interpolated_dihedral(ls_dihedrals : np.ndarray, dihr_of_interest : float) -> float :
    """
    This function is well important!
    Receives the dihedrals calculated when using dihedral_array(), where we calculate the angle
        between three atoms and the cone of vectors, which represent the atom of interest.
    Then we use the value we know is the correct dihedral angle and we check which in between
        which vectors the value of interest is.
    We then use the index of the two vectors and retrieve the theta angle, because we rotated the cone
        in a certain way, but the order still matters. Theta is used to generate the cone vector
        With the two theta values in place, we interpolate between the theta and the calculated dihedral angles.
        The returned value theta is the exact theta angle we need to use in order to generate
            the single_vector that already has the correct phi angle.
    """

    # Check the values of the range of dihedrals and the value of interest
    #print(ls_dihedrals) ; print(dihr_of_interest)

    # Generate an array of theta values like when we generate them
    theta = np.linspace(0, 2*np.pi, num=18, endpoint=False)

    # Determine if val(ls_dihedrals) is ascending or descending
    slope = LFT2.check_slope_of_array(ls_dihedrals)

    if slope == "DESCENDING":
        ## Let's check the difference in between ls_dihedrals[i] and 180 or -180 is less than 20
        for i in range(len(ls_dihedrals)):

            # if the value[i] is close to 180, check if the value is less than 20 
            if 180.0 - ls_dihedrals[i] <= 20.0:
                if ls_dihedrals[i] <= dihr_of_interest <= 180.0:
                    lowerRAD = theta[i]

                    # upperbound is difference between i en 180, converted to radians and added to lowerbound
                    upperRAD = theta[i] - ((180.0 - ls_dihedrals[i]) * (np.pi/ 180.0))
                    dihr_boundaries = ((lowerRAD, ls_dihedrals[i]), (upperRAD, 180.0))
                    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

            # if the value[i] is close to -180, the difference is less than -340 (-180 + -180 = -360)
            if -180.0 + ls_dihedrals[i] <= -340.0:
                if  -180.0 <= dihr_of_interest <= ls_dihedrals[i]:
                    lowerRAD = theta[i]

                    # upperbound is difference between i and -180, convertedto radians and added to lowerbound
                    upperRAD = theta[i] - ((-180 - ls_dihedrals[i]) * (np.pi/180.0))
                    dihr_boundaries = ((lowerRAD, ls_dihedrals[i]), (upperRAD, -180.0))
                    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

            # if it's the last index, compare between last and first index, to make it circular
            if i == (len(ls_dihedrals) - 1):
                if ls_dihedrals[i] <= dihr_of_interest <= ls_dihedrals[0]:
                    dihr_boundaries = ((theta[i], ls_dihedrals[i]), (theta[0], ls_dihedrals[0]))

                    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

            # if it's not last index, just carry on as usual
            if ls_dihedrals[i] >= dihr_of_interest >= ls_dihedrals[i+1]:
                dihr_boundaries = ((theta[i], ls_dihedrals[i]), (theta[i+1], ls_dihedrals[i+1]))
                #print(dihr_boundaries)

                return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

    if slope == "ASCENDING":
        for i in range(len(ls_dihedrals)):

            # if the value[i] is close to 180, check if the value is then between value[i] and 
            if 180 - ls_dihedrals[i] <= 20:
                if ls_dihedrals[i] <= dihr_of_interest <= 180:
                    lowerRAD = theta[i]

                    # upperbound is diff between i and 180, converted to radians and added to lowerbound
                    upperRAD = theta[i] - ((180.0 - ls_dihedrals[i]) * (np.pi/ 180.0))
                    dihr_boundaries = ((lowerRAD, ls_dihedrals[i]), (upperRAD, 180.0))
                    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

            # if the value[i] is close to -180, the difference is less than -340 (-180 + -180 = -360)
            if -180.0 + ls_dihedrals[i] <= -340.0:
                if  -180.0 <= dihr_of_interest <= ls_dihedrals[i]:
                    lowerRAD = theta[i]

                    # upperbound is difference between i and -180, convertedto radians and added to lowerbound
                    upperRAD = theta[i] - ((-180 - ls_dihedrals[i]) * (np.pi/180.0))
                    dihr_boundaries = ((lowerRAD, ls_dihedrals[i]), (upperRAD, -180.0))
                    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

            # if it's the last index, compare between last and first index, to make it circular
            # NB: THIS DIFFERS FROM THE FOR LOOP ABOVE, WE GO IN THE OPPOSITE DIRECTION
            if i == (len(ls_dihedrals) - 1):
                if ls_dihedrals[0] <= dihr_of_interest <= ls_dihedrals[i]:
                    dihr_boundaries = ((theta[i], ls_dihedrals[i]), (theta[0], ls_dihedrals[0]))

                    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

            # if it's not last index, just carry on as usual
            # NB: THIS DIFFERS FROM THE FOR LOOP ABOVE, WE GO IN THE OPPOSITE DIRECTION
            if ls_dihedrals[i] <= dihr_of_interest <= ls_dihedrals[i+1]:
                dihr_boundaries = ((theta[i], ls_dihedrals[i]), (theta[i+1], ls_dihedrals[i+1]))
                #print(dihr_boundaries)

                return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)


def position_linker(O5_nucleoAtom : np.array, P_vector : np.array, linker : np.ndarray) -> np.ndarray:
    """ Get the translation vector of the linker and then move.
    Returns the position at which we move the linker segment to.
    """
    # Stretch the normalized vector; get the position of the P_atom
    p_position = O5_nucleoAtom + (P_vector * 1.6)

    #Get the distance required to move the linker to the P_atom position
    p_move_to = p_position - linker.get_P()

    # move the entire linker to the P_atom position
    return p_move_to + linker.array


def move_vector_to_loc(array_of_molecules : np.ndarray, distance_to_loc : np.array) -> np.ndarray:
    """ Move the array of vectors from the origin to the place you want to have it"""
    return array_of_molecules + distance_to_loc


def move_vector_to_origin(array_of_molecules : np.array, distance_to_origin : np.array) -> np.ndarray:
    """ Move the vector to the origin by the distance of a specific atom to that origin """
    return array_of_molecules - distance_to_origin


def create_PDB_from_matrix(matrix : np.ndarray, nucleoside : np.ndarray, linker : np.ndarray) -> None:
    print("Writing to pdb ...")

    df = pd.DataFrame()

    df['RecName'] = ['ATOM' for x in range(matrix.shape[0])]
    df['AtomNum'] = np.arange(start=1, stop=matrix.shape[0] + 1)
    df['AtomName'] = json.loads(nucleoside.jason['pdb_properties']['Atoms']) + json.loads(linker.jason['pdb_properties']['Atoms']) + json.loads(nucleoside.jason['pdb_properties']['Atoms'])
    df['AltLoc'] = ' '
    df['ResName'] = 'POS'     # This is temporary, just to see what the results are
    df['Chain'] = 'A'
    df['Sequence'] = str(1)
    df['X_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,0]))
    df['Y_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,1]))
    df['Z_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,2]))
    df['Occupancy'] = '1.00'
    df['Temp'] = '0.00'
    df['SegmentID'] = str('   ')
    df['ElementSymbol'] = json.loads(nucleoside.jason['pdb_properties']['Symbol']) + ['P', 'O', 'O'] + json.loads(nucleoside.jason['pdb_properties']['Symbol'])

    with open('testing_daedalus.pdb' ,'w') as pdb:
        for index, row in df.iterrows():
            split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
            pdb.write('%-6s%5s%5s%s%s%2s%5s   %8s%8s%8s%6s%6s%4s      %2s\n' % tuple(split_line))
        pdb.write('END')

##---------------------------- EULER ANGLE ROTATION MATRIX, NOT USED ANYMORE ----------------##
#def generate_and_rotate_single_vector_EULER(interpolated_theta_angle, phi, rotation_matrix):
#
#    ### PHI NEEDS TO BE IN RADIANS
#    single_vector = np.array([ 1.0 * np.cos(interpolated_theta_angle) * np.sin(phi),
#                               1.0 * np.sin(interpolated_theta_angle) * np.sin(phi),
#                               1.0 * np.cos(phi)]
#                            ).T
#
#    rotated_vector = np.dot(rotation_matrix, single_vector)
#
#    return rotated_vector
#
#def get_rM(vector_to_rotate_onto, vector_to_rotate_from=np.array([0,0,1])):
#    """
#    EULER ANGLES
#    Find the rotation matrix (rM) associated with counterclockwise rotation
#    about the given axis by theta radians.
#    Credit: http://stackoverflow.com/users/190597/unutbu
#
#    Args:
#        axis (list): rotation axis of the form [x, y, z]
#        theta (float): rotational angle in radians
#
#    Returns:
#        array. Rotation matrix.
#    """
#
#    # To rotate a matrix, we make a dot product with the rotation matrix and the vector we want to move.
#    #       The results is the rotated matrix. 
#
#    # The cone's axis is vector_to_rotate_from
#    # np.array([0, 0, 1])
#
#    axis = get_direction_for_rM(vector_to_rotate_from, vector_to_rotate_onto)      # DIRECTION
#    theta = get_angle_for_rM(vector_to_rotate_from, vector_to_rotate_onto)         # ANGLE
#
#    axis = axis/np.sqrt(np.dot(axis, axis))
#    a = np.cos(theta/2.0)
#    b, c, d = -axis*np.sin(theta/2.0)
#    aa, bb, cc, dd = a*a, b*b, c*c, d*d
#    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
#    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
#                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
#                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
#
#
#def rotate_cone_vector(rotation_matrix, cone_vector):
#    # Create empty vector
#    rotated_cone = np.empty(shape=(cone_vector.shape[0], cone_vector.shape[1]))
#    for i in range(cone_vector.shape[0]):
#        rotated_cone[i] = np.dot(rotation_matrix, cone_vector[i])
#
#    return rotated_cone
#
#
#def rotate_single_vector(rotation_matrix, atoms):
#    # rotate a single 
#    return np.dot(rotation_matrix, atoms)