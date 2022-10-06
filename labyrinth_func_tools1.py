import numpy as np
import numpy.linalg as LA
from scipy.spatial.transform import Rotation as R

from typing import Tuple
import json
import labyrinth_func_tools2 as LFT2

"""     labyrinth_func_tools1 contains all the mathematical functions required, in python, to calculate for the rotations.  """

def return_normalized(vector : np.ndarray) -> np.ndarray:
    """ returns a normalized vector """
    return vector / LA.norm(vector)


def generate_cone_vector(phi_angle : float) -> np.ndarray:
    """ Phi is the angle the vector makes with the Z-axis
        Graphing Spherical Coordinates in GeoGebra 3D (Part 2): A Cone about z-axis
            https://www.youtube.com/watch?v=Wl3Z3AqfI6c

        Generating uniform unit random vectors in R^n
        Prompted phi in rad.
        conical vector = array[ones] * array[theta] * array[phi]    """

    theta = np.linspace(0, 2*np.pi, num=18, endpoint=False)
    phi = np.full(18, phi_angle)

    return np.array([ np.cos(theta) * np.sin(phi),
                      np.sin(theta) * np.sin(phi),
                      np.cos(phi)]
                    ).T


def get_direction_of_rotation(from_vector : np.ndarray, vector_to_rotate_onto : np.ndarray) -> np.ndarray:
    """ Cross product with the cone's axis to get the direction of the rotation axis
        Get the direction where vectorA rotates onto vectorB ; u = vectora X vectorb
        Let's normalize the direction """

    cross_product = np.cross(from_vector, vector_to_rotate_onto)
    return return_normalized(cross_product)


def get_angle_of_rotation(from_vector : np.ndarray, vector_to_rotate_onto : np.ndarray) -> float:
    """ The scalar product (dot product) to get the cosine angle.
        Here we do the arccos, to get the angle immediately. """
    return np.arccos(np.dot(from_vector, vector_to_rotate_onto))


def get_length_of_vector(vector1 : np.array, vector2 : np.array) -> float:
    """ Uses the linear algebra module to return the length of the vector """
    return LA.norm(vector2 - vector1)


def TESTTHETAget_quaternion(vector_to_rotate_onto : np.ndarray, vector_to_rotate_from : np.array = np.array([0,0,1])) -> np.array:
    """ Quaternion mathematics using the scipy.spatial.transform.Rotation library
        Create the quaternion that is associated with the angle and axis of rotation
        Apply it to the vector you want to rotate.  """
    # axis has been normalised
    axis = get_direction_of_rotation(vector_to_rotate_from, vector_to_rotate_onto)      # DIRECTION

    # angle already in radians 
    theta = np.pi / 4     # ANGLE

    # Create the quaternion
    qx = axis[0] * np.sin(theta)
    qy = axis[1] * np.sin(theta)
    qz = axis[2] * np.sin(theta)
    qw = np.cos(theta)

    quaternion = R.from_quat([qx, qy, qz, qw])

    return quaternion


def get_quaternion(vector_to_rotate_onto : np.ndarray, vector_to_rotate_from : np.array = np.array([0,0,1])) -> np.array:
    """ Quaternion mathematics using the scipy.spatial.transform.Rotation library
        Create the quaternion that is associated with the angle and axis of rotation
        Apply it to the vector you want to rotate.  """
    # axis has been normalised
    axis = get_direction_of_rotation(vector_to_rotate_from, vector_to_rotate_onto)      # DIRECTION

    # angle already in radians 
    theta = get_angle_of_rotation(vector_to_rotate_from, vector_to_rotate_onto) / 2     # ANGLE

    # Create the quaternion
    qx = axis[0] * np.sin(theta)
    qy = axis[1] * np.sin(theta)
    qz = axis[2] * np.sin(theta)
    qw = np.cos(theta)

    quaternion = R.from_quat([qx, qy, qz, qw])

    return quaternion

def get_quaternion_custom_axis(axis : np.array, theta : float) -> np.array:
    """ Generate a customized quaternion. This is only succesful if the direction axis is perpendicular. """
    # Create the quaternion
    qx = axis[0] * np.sin(theta)
    qy = axis[1] * np.sin(theta)
    qz = axis[2] * np.sin(theta)
    qw = np.cos(theta)

    quaternion = R.from_quat([qx, qy, qz, qw])

    return quaternion


def get_normal_vector_of_plane(v_first : np.ndarray, v_second : np.ndarray) -> np.ndarray:
    """ Retrieve the normal of a plane, by taking the cross product of two vectors. """
    v_first = return_normalized(v_first)
    v_second = return_normalized(v_second)
    cross_product = np.cross(v_first, v_second)

    return cross_product / LA.norm(cross_product)


def rotate_with_quaternion(quaternion, vector : np.ndarray) -> np.ndarray:
    """ Vector rotation through quaternion mathematics"""

    return quaternion.apply(vector)


def generate_and_rotate_single_vector(interpolated_theta_angle : float, phi : float, quaternion) -> np.ndarray:
    """ Generate a single vector with the correct angle.
        Then rotate said angle to the correct orientation using the previously used quaternion.
        As always, phi needs to be in RADIANS
        vector = r* theta * phi     """

    single_vector = np.array([ np.cos(interpolated_theta_angle) * np.sin(phi),
                               np.sin(interpolated_theta_angle) * np.sin(phi),
                               np.cos(phi)]
                            ).T

    rotated_vector = quaternion.apply(single_vector)

    return rotated_vector


def assert_dot_product(_dotProduct : np.array) -> bool :
    """ We have noticed that whenever two vectors (i.e. j1 and j2, where j1 is imposed upon j2) make a dot product of
        -1 or -180 degrees, that the rotation of the quaternion spazzes out. To combat this, we will test whether or not
        j1 and j2 have a dot product nearing this value.

        If this is the case, we simply alter the orientation of the initial array by 45 degrees (arbritrary value) and then return this and reposition
        the nucleoside array.

        if the difference(-1, _dotProduct) is close to zero, then change the orientation. """
    _testAgainst = - 1 - _dotProduct
    #print(_testAgainst)

    # Check if the asserted dot product is close to zero, from the negative range up
    if -0.1 <= _testAgainst <= 0 :
        return True
    return False


def reorient_nucleoside_array(complementary_array : np.array) -> np.ndarray :
    """ Reorient the array that has the weird, gimmicky dotproduct equals -1 thing from up above.

        What we did here is choose the standard x axis to rotate over, bring the nucleoside array by a random atom to the origin
        and rotate it, only to put it back to its original location."""

    # Custom axis and custom angle to rotate. We rotate over an angle of 45 degrees just because it rotates it enough to steer away from the point
    _quaternion_custom = get_quaternion_custom_axis(np.array([1,0,0]), (np.pi/4))

    return move_to_origin_ROTATE_move_back_to_loc(_quaternion_custom, complementary_array, complementary_array[0])

#def apply_rotation_of_planes(quaternion : np.array, nucleoside_array : np.ndarray, index_of_vectors : list) -> np.ndarray :
#    """ This functions applies a second rotation to the nucleoside. The problem here is that we cannot use a custom axis, since direction axises are always perpendicular to the two vectors.
#        This is why we rotate the array as usual. But we only override the vectors that we actually want to rotate.
#        This way, we effectively rotate over a certain bond, but the mathematical part just did a copy-paste of the vectors we really want
#        index_of_vectors : the vectors (or the atoms) that we do not want to have rotated. These will be cut out the list of indices.
#        """
#    # Get size of the array, then make a list of the range of indices to emulate the indices of the array we calculate with
#   # size_of_array = nucleoside_array.shape[0]
#   # list_of_indices = [i for i in range(size_of_array)]
#
#    # Remove index of the vectors of array that we do not want to rotate, so we do not override them later in this function
#    #for index in index_of_vectors:
#    #    list_of_indices.remove(index)
#
#    # apply rotation
#    #tmp_rotated_vector = quaternion.apply(nucleoside_array)
#    return quaternion.apply(nucleoside_array)
#
#    # Override the vectors into the nucleoside_array, ofcourse already having ommited the vectors you do not want to override
#   # for i in list_of_indices:
#   #     nucleoside_array[i] = tmp_rotated_vector[i]
#
#   # return nucleoside_array
#

def dihedral_array(atoms_in_dihr : np.ndarray, cone_vector : np.ndarray) -> np.ndarray :
    """ https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    Praxeolitic dihedral calculation. as of now, the most efficient way to calculate for a dihedral angle.  """
    # Example dihedral : C4' - C5' - O5' - P

    v2 = atoms_in_dihr[0] # O5'
    v1 = atoms_in_dihr[1] # C5'
    v0 = atoms_in_dihr[2] # C4'
    v3 = cone_vector    # P , this contains all the possible P atoms

    b0 = -1.0 * (v1 - v0)   #C4-C5
    b1 = v2 - v1            #C5-O5
    b2 = v3                 #O5-P

    # normalize b1 so that it does not influence magnitude of vector rejections that come next
    b1 = return_normalized(b1)

    ## vector rejections
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


def dihedral_single(first : float, second : float, third : float, fourth : float) -> float:
    """ https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    Praxeolitic dihedral calculation. as of now, the most efficient way to calculate for a dihedral angle.  """
    # Example dihedral : C4' - C5' - O5' - P

    v2 = third              # O5'
    v1 = second             # C5'
    v0 = first              # C4'
    v3 = fourth             # Supposed P

    b0 = (v1 - v0) * -1.0   # C4' - C5' 
    b1 = v2 - v1            # C5' - O5'
    b2 = v3 - v2            # O5' - P

    # normalize b1 so that it does not influence magnitude of vector rejections that come next
    b1 = return_normalized(b1)

    ## vector rejections
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return np.degrees(np.arctan2(y, x))


def interpolate_dihedrals(tuple_dihr : Tuple[Tuple[float, float], Tuple[float, float]], angle_dihr : float) -> float:
    """     y = y1 + [ (x - x1) / (x2 - x1) * (y2 - y1) ]

    tuple_dihr = ( (rad1, degrees1), (rad2, degrees2) )

    The y's are the theta angles, indexed by the keys of the dictionary being parsed.    """
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
                if ls_dihedrals[i] >= dihr_of_interest >= ls_dihedrals[0]:
                    dihr_boundaries = ((theta[i], ls_dihedrals[i]), (2 * np.pi, ls_dihedrals[0]))

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
                    dihr_boundaries = ((theta[i], ls_dihedrals[i]), (2 * np.pi, ls_dihedrals[0]))

                    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)

            # if it's not last index, just carry on as usual
            # NB: THIS DIFFERS FROM THE FOR LOOP ABOVE, WE GO IN THE OPPOSITE DIRECTION
            if ls_dihedrals[i] <= dihr_of_interest <= ls_dihedrals[i+1]:
                dihr_boundaries = ((theta[i], ls_dihedrals[i]), (theta[i+1], ls_dihedrals[i+1]))
                #print(dihr_boundaries)

                return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)


def move_vector_to_loc(array_of_molecules : np.ndarray, distance_to_loc : np.array) -> np.ndarray:
    """ Move the array of vectors from the origin to the place you want to have it"""
    return array_of_molecules + distance_to_loc


def move_vector_to_origin(array_of_molecules : np.array, distance_to_origin : np.array) -> np.ndarray:
    """ Move the vector to the origin by the distance of a specific atom to that origin """
    return array_of_molecules - distance_to_origin


def move_to_origin_ROTATE_move_back_to_loc(quaternion : np.ndarray, array_to_manipulate : np.ndarray, distance_to_move : np.array) -> np.ndarray:
    """ Takes out the iterative process of moving an array to the origin, rotating it and moving it back again to it's original spot """
    # move to origin
    array_to_origin = move_vector_to_origin(array_to_manipulate, distance_to_move)
    # rotate
    rotated_array = rotate_with_quaternion(quaternion, array_to_origin)
    # move back to the original location
    return move_vector_to_loc(rotated_array, distance_to_move)


def assert_length_of_vector(length : float) -> bool:
    """ Check if length is roughly the correct size or if the nucleotide needs to be rotated """
    return 1.30 <= length <= 1.90


def assert_size_of_angle(angle : float, angle_to_fit : float) -> bool:
    """ Check if angle is roughly the correct size or if the nucleotide needs to be rotated """
    return angle_to_fit * 0.8 <= angle <= angle_to_fit * 1.2


def assert_dihedral_of_nucleotide(calculated_dihedral : float, dihedral : float) -> bool:
    """ Check if dihedral is roughly the correct angle or if the nucleotide needs to be rotated """
    offset_dihedral = 25
    # If the dihedral surpasses passed 180 degrees, then assert if the dihedral is smaller than -180 + the difference of (180 - dihedral)
    if dihedral + offset_dihedral > 180:
        diff = -180 + (dihedral - 180 + offset_dihedral)

        check1 = dihedral - offset_dihedral <= calculated_dihedral <= 180
        check2 = -180 <= calculated_dihedral <= diff

        if check1 or check2:
            return True

    # If the dihedral goes below -180 degrees, then assert if the dihedral is smaller than 180 - the difference of (-180 + dihedral)
    if dihedral - offset_dihedral < -180:
        diff = 180 + (calculated_dihedral + 180 - offset_dihedral)

        check1 = -180 <= calculated_dihedral <= dihedral + offset_dihedral
        check2 = diff <= calculated_dihedral <= 180

        if check1 or check2:
            return True

    # If it does not go over the bounds of -180 or 180, then just calculate it like this
    return dihedral - offset_dihedral <= calculated_dihedral <= dihedral + offset_dihedral


def smallest_difference(array : np.array, ref_value : float) -> np.array:
    """ Returns the value in the array where the difference with the reference value is the smallest"""
    diff_array = array - ref_value
    return np.abs(diff_array)

#def check_phi_angle_of_vector(vectors : np.ndarray, axis : np.array = np.array([0,0,1])) :
#    """ The dotproduct determines the angle of the vector with a given axis/vector.
#    This function is mainly for debugging purposes and has no value for building duplexes.
#
#    If the angle is not what you expected it would be, try inverting the direction of one of the vectors.
#
#    EDIT : This process has been automated and this function is redudant, but I found it useful when I learned about linear algebra
#        to see where I was going wrong so it STAYS!"""
#
#    # if a vector is given and it is not k^-hat
#    if list(axis) != [0,0,1]:
#        arr = np.empty(shape=(vectors.shape[0]))
#
#        for i in range(vectors.shape[0]):
#            arr[i] = np.arccos(np.dot(vectors[i], axis))
#
#        #print(np.degrees(arr))
#        arr = np.around(arr, decimals=10)                   # Round of values at the end to avoid errors on equality at the 17th decimal place
#        if list(arr).count(arr[0]) != len(arr):
#            print('Not all cone angles are correctly generated')
#            sys.exit(0)
#
#        return
#
#    # If no vector is given, default to k^-hat; Z_axis)
#    arr = np.arccos(np.dot(vectors, axis))
#    #print(np.degrees(arr))
#    arr = np.around(arr, decimals=10)
#    if not list(arr).count(arr[0]) == len(arr):
#
#        print('Not all cone angles are correctly generated')
#        sys.exit(0)
#
