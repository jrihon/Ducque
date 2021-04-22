"""

"""
import numpy as np
import json
import sys

# CLASSES
class Nucleoside:

    def __init__(self,jsonfile):
        self.splitted = jsonfile.split('.')[0]

        with open(jsonfile, 'r') as jsonf:
            self.jason = json.load(jsonf)

        self.array =  np.asarray(json.loads(self.jason['pdb_properties']['Coordinates']), dtype=float)


    def normalize_vector(self, vector):
        # Retrieve the magnitude of the vector when normalized.
        # Return the vector that has been normalized, so ... rework this piece of code

        normVec = np.linalg.norm(self.array[0])
        return vector / normVec

    def get_beta(self):
        # because beta is still inside a the backbone dictionary, we need to load it again 
        return float(json.loads(self.jason['Dihedrals']['Backbone'])['beta'])


class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_COP(self):
        # since the angles are inside the angles dictionary, we don't need to load string again
        return float(self.jason['Angles']['C5_O5_P'])

# FUNCTIONS
def generate_cone_vector(phi_angle):
    " Phi is the angle the vector makes with the Z-axis """

    """
    Graphing Spherical Coordinates in GeoGebra 3D (Part 2): A Cone about z-axis
    https://www.youtube.com/watch?v=Wl3Z3AqfI6c

    Generating uniform unit random vectors in R^n
    """
    ### PHI_ANGLE NEEDS TO BE IN RADIANS
    rho = np.full(18, 1)
    theta = np.linspace(0, 2*np.pi, num=18, endpoint=False)
    phi = np.full(18, phi_angle)

    return np.array([ rho * np.cos(theta) * np.sin(phi),
                      rho * np.sin(theta) * np.sin(phi),
                      rho * np.cos(phi)]
                    ).T

def generate_rotate_single_vector(rotation_matrix, interpolated_theta_angle, phi):

    ### PHI NEEDS TO BE IN RADIANS
    single_vector = np.array([ 1 * np.cos(interpolated_theta_angle) * np.sin(phi),
                              1 * np.sin(interpolated_theta_angle) * np.sin(phi),
                              1 * np.cos(phi)]
                            ).T

    rotated_vector = np.dot(rotation_matrix, single_vector)

    return rotated_vector


def get_direction_for_rM(from_vector, vector_to_rotate_onto):

    # cross product with the cone's axis to get the direction of the rotation axis
    # Get the direction where vectorA rotates onto vectorB ; u = vectora X vectorb
    #       Let's normalize the direction

    return np.cross(from_vector, vector_to_rotate_onto) / np.linalg.norm(np.cross(from_vector, vector_to_rotate_onto))

def get_angle_for_rM(from_vector, vector_to_rotate_onto):

    # The scalar product (dot product) to get the cosine angle. Here we do the arccos, the get the angle immediately.

    return np.arccos(np.dot(from_vector, vector_to_rotate_onto))

def get_rM(vector_to_rotate_onto):
    """
    Find the rotation matrix (rM) associated with counterclockwise rotation
    about the given axis by theta radians.
    Credit: http://stackoverflow.com/users/190597/unutbu

    Args:
        axis (list): rotation axis of the form [x, y, z]
        theta (float): rotational angle in radians

    Returns:
        array. Rotation matrix.
    """

    # To rotate a matrix, we make a dot product with the rotation matrix and the vector we want to move.
    #       The results is the rotated matrix. 

    # The cone's axis
    _m = np.array([0, 0, 1])

    axis = get_direction_for_rM(_m, vector_to_rotate_onto)      # DIRECTION
    theta = get_angle_for_rM(_m, vector_to_rotate_onto)         # ANGLE

    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def rotate_cone_vector(rotation_matrix, cone_vector):
    # Create empty vector
    rotated_cone = np.empty(shape=(cone_vector.shape[0], cone_vector.shape[1]))
    for i in range(cone_vector.shape[0]):
        rotated_cone[i] = np.dot(rotation_matrix, cone_vector[i])

    return rotated_cone


def check_phi_angle_of_vector(vectors, axis = np.array([0,0,1])):
    """ The dotproduct determines the angle of the vector with a given axis/vector """

    # if a vector is given and it is not k^-hat
    if list(axis) != [0,0,1]:
        arr = np.empty(shape=(vectors.shape[0]))

        for i in range(vectors.shape[0]):
            arr[i] = np.arccos(np.dot(vectors[i], axis))

        if not list(arr).count(arr[0]) == len(arr):
            print('Not all cone angles are correctly generated')
            sys.exit(0)

        return

    # If no vector is given, default to k^-hat; Z_axis)
    arr = np.degrees(np.arccos(np.dot(vectors, axis)))
    if not list(arr).count(arr[0]) == len(arr):
        print('Not all cone angles are correctly generated')
        sys.exit(0)


def praxeolitic_dihedralRANGE(json_array, cone_vector):
 
    """
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    """

    v2 = json_array[0] # O5'
    v1 = json_array[1] # C5'
    v0 = json_array[4] # C4'
    v3 = cone_vector    # P , this contains all the possible P atoms

    b0 = -1.0 * (v1 - v0)   #C4-C5
    b1 = v2 - v1            #C5-O5
    b2 = v3                 #O5-P

    # normalize b1 so that it does not influence magnitude of vector rejections that come next
    b1 /= np.linalg.norm(b1)

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


def praxeolitic_dihedralSINGLE(json_array, single_vector):

    v2 = json_array[0] # O5'
    v1 = json_array[1] # C5'
    v0 = json_array[4] # C4'
    v3 = single_vector # Supposed P

    b0 = -1.0 * (v1 - v0)
    b1 = v2 - v1 #C5O5
    b2 = v3

    # normalize b1 so that it does not influence magnitude of vector rejections that come next
    b1 /= np.linalg.norm(b1)

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


def interpolate_dihedrals(dict_dihr, angle_dihr):
    """
    y = y1 + [ (x - x1) / (x2 - x1) * (y2 - y1) ]
    """

    theta = np.linspace(0, 2*np.pi, num=18, endpoint=False)
    # the y's are the theta angles, indexed by the keys of the dictionary being parsed 
    y2 = theta[ int(list(dict_dihr.keys())[1]) ]
    y1 = theta[ int(list(dict_dihr.keys())[0]) ]
    x2 = list(dict_dihr.values())[1]
    x1 = list(dict_dihr.values())[0]
    x = angle_dihr

    return y1 + ( (x - x1) / (x2 - x1) * (y2 - y1) )


def get_interpolated_dihedral(ls_dihedrals, dihr_of_interest):

    # Create dictionary to assume values of boundaries and the calculated angle
    for i in range(len(ls_dihedrals)):
        if i == (len(ls_dihedrals) - 1):
            if ls_dihedrals[i] < dihr_of_interest < ls_dihedrals[0]:
                dihr_boundaries = { str(i) : ls_dihedrals[i], str(0) : ls_dihedrals[0] }
                break
        if ls_dihedrals[i] < dihr_of_interest < ls_dihedrals[i+1]:
            dihr_boundaries = { str(i) : ls_dihedrals[i], str(i+1) : ls_dihedrals[i+1] }
            break

    # The interpolate function exists in this module
    return interpolate_dihedrals(dihr_boundaries, dihr_of_interest)



