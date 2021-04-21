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


def get_rotation_matrix(axis, theta):
    """
    Find the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
    Credit: http://stackoverflow.com/users/190597/unutbu

    Args:
        axis (list): rotation axis of the form [x, y, z]
        theta (float): rotational angle in radians

    Returns:
        array. Rotation matrix.
    """

    #axis = np.asarray(axis)
    #theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


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


def Architecture(nucleic_acid):

    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
    nucleo = labyrinth_func.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

    # Calculate the normal of C4'-C5'-O5' 
    v0 = nucleo.array[0] # O5'
    v1 = nucleo.array[1] # C5'
    v2 = nucleo.array[4] # C4'

    C5O5 = (v1 - v0) * -1.0        # I-hat , where origin is equal to v1
    C5C4 = v2 - v1                  # J-hat

    # Get vector of normalized magnitude 
    C5O5 = C5O5 / np.linalg.norm(C5O5)

    # angle between C5' - O5' - P in radians, convert from degree to radian
    _angleCOP = 118.980 * (np.pi / 180)

    """
    Graphing Spherical Coordinates in GeoGebra 3D (Part 2): A Cone about z-axis
    https://www.youtube.com/watch?v=Wl3Z3AqfI6c

    Generating uniform unit random vectors in R^n
    """

    rho = np.full(18, 1)
    theta = np.linspace(0, 2*np.pi, num=18, endpoint=False)
    phi = np.full(18, _angleCOP)

    vectors_cone = np.array([ rho * np.cos(theta) * np.sin(phi),
                                        rho * np.sin(theta) * np.sin(phi),
                                        rho * np.cos(phi)]
                                         ).T

    # Check if the dotproduct is the same everywhere
    _mcheck = np.array([0, 0, 1])
    _rotcheck = np.arccos(np.dot(vectors_cone, _mcheck))    # this gives the correct angles to the vectors


    # Rotate the vector fit over the C505 vector we want to discuss
    # find rotation axis _u and toration angle _rot. We rotate around the z-axis here ' the cone axis '.

    # To rotate a matrix, we make a dot product with the rotation matrix and the vector we want to move.
    #       The results is the rotated matrix. 

    # The cone's axis
    _m = np.array([0, 0, 1])

    # cross product with the cone's axis to get the direction of the rotation axis
    # Get the direction where vectorA rotates onto vectorB ; u = vectora X vectorb
    #       Let's normalize the direction
    _uDirection = np.cross(_m, C5O5) / np.linalg.norm(np.cross(_m, C5O5))

    # The scalar product (dot product) to get the cosine angle. Here we do the arccos, the get the angle immediately.
    _rotAngle = np.arccos(np.dot(_m, C5O5))

    # Get the rotation matrix, which will be multiplied with the vector
    _rot_matrix = get_rotation_matrix(_uDirection, _rotAngle)

    rotated_cone = np.empty(shape=(vectors_cone.shape[0], vectors_cone.shape[1]))
    for i in range(vectors_cone.shape[0]):
        rotated_cone[i] = np.dot(_rot_matrix, vectors_cone[i])

    #------------------------------------------------------------------------------------ Cone is in the correct location
    # Prints the angle between the cone vectors and the C505 vector. Only required to check if it's correct
    dot_vector = np.empty(shape=(vectors_cone.shape[0]))
    for i in range(vectors_cone.shape[0]):
        dot_vector[i] = np.dot(C5O5, rotated_cone[i])

    #--------------------------------------------------------------- Praxeolitic formula
    """
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    
    """

    v2 = nucleo.array[0] # O5'
    v1 = nucleo.array[1] # C5'
    v0 = nucleo.array[4] # C4'
    v3 = rotated_cone    # P , this contains all the possible P atoms

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
    dihedrals = np.empty(b2.shape[0])

    for i in range(b2.shape[0]):
        v = b0 - np.dot(b0, b1)*b1
        w = b2[i] - np.dot(b2[i], b1)*b1

        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        dihedrals[i] = np.degrees(np.arctan2(y, x))
#------------------------------------------------------------ Now that we have a variety of dihedrals, we try to interpolate
    # DIHEDRAL ANGLE TO FIT TO
    beta_dihr = float(json.loads(nucleo.jason['Dihedrals']['Backbone'])['beta'])

    for i in range(len(dihedrals)):
        if i == (len(dihedrals) - 1):
            if dihedrals[i] < beta_dihr < dihedrals[0]:
                angle_extremities = { str(i) : dihedrals[i], str(0) : dihedrals[0] }
                break
        if dihedrals[i] < beta_dihr < dihedrals[i+1]:
            angle_extremities = { str(i) : dihedrals[i], str(i+1) : dihedrals[i+1] }
            break

    y_interpolate = interpolate_dihedrals(angle_extremities, beta_dihr)

#-------------------------------------------------------------------- use interpolated value to check if the dihedral is correct
    checkVector = np.array([ 1 * np.cos(y_interpolate) * np.sin(_angleCOP),
                                        1 * np.sin(y_interpolate) * np.sin(_angleCOP),
                                        1 * np.cos(_angleCOP)]
                                         ).T

    rotated_vector = np.dot(_rot_matrix, checkVector)


    #print(np.degrees(np.dot(C5O5, rotated_vector)))

    v2 = nucleo.array[0] # O5'
    v1 = nucleo.array[1] # C5'
    v0 = nucleo.array[4] # C4'
    v3 = rotated_vector  # P , this contains all the possible P atoms

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
    np.degrees(np.arctan2(y, x))
    #print(beta_dihr)
    #print(np.degrees(np.arctan2(y, x)))
    # Parse dictionary for the correct linker segment
    #link = labyrinth_func.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#
#    # vector to rotate the z_axis around
#    #ax.scatter(C5O5[0], C5O5[1], C5O5[2], color='blue', alpha=0.5)
#    ax.scatter(b0[0], b0[1], b0[2], color='blue')
#    ax.scatter(b1[0], b1[1], b1[2], color='orange')
#    ax.scatter(b2[:,0], b2[:,1], b2[:,2], color='red')
#    ax.scatter(b2[16,0], b2[16,1], b2[16,2], color='green')
#    #rotated cone
#    #ax.scatter(rotated_cone[:,0], rotated_cone[:,1], rotated_cone[:,2], color='green', alpha=0.5)
#
#    #original cone
#    #ax.scatter(vectors_cone[:,0], vectors_cone[:,1], vectors_cone[:,2], color='r', alpha=0.5)
#
#    #origin
#    ax.scatter(0,0,0, color='black')
#
#    ax.set_zlim(-2,2)
#    ax.set_xlim(-2,2)
#    ax.set_ylim(-2,2)
#    ax.set_xlabel('x_axis')
#    ax.set_ylabel('y_axis')
#    ax.set_zlabel('z_axis')
#
#    plt.show()
