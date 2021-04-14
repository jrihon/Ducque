"""
Godot game engine that simplifies how translation and rotation works.
https://docs.godotengine.org/en/2.1/learning/features/math/vector_math.html

http://www.sciencebits.com/dot_product

Vector Math is, in a great deal, dimension-independent so adding or removing an axis only adds very little complexity!

Some formulas
Consider a = np.array([5,3,7])
Consdier h = np.array([6,9,3])

Length of a vector:
    length_a = np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])

Angle from a vector:
    angle_a = np.atan2(a)
    This is used to calculate the angle of the vector with... an axis that has to be predefined?

Rotate a 2D vector:
    swap the x and y, then negate either (depending on the turn) coordinates for a 90 degree turn. Or 180 degree turn and negate both.

Normalise a vector:
    norm_a = np.linalg.norm() ; in which you divide the coordinates of the vector by the length of the vector. 

Dot product:
    the dot product returns a scalar.
    dot_prod = a[0]*h[0] + a[1]*h[1] + a[2]*h[2]

    The order of addition does not matter so that eventually so that np.dot(a,b) and np.dot(b,a) is fine. Not to confuse with matrix multiplication, which is completely different.

    Remember; a scalar is just a single floating number (float == scalar)

    ANGLE
    If dot_prod > 0 (greater than), that means that the angle between the two vectors is less than (<) 90 degrees
    If dot_prod == 0 (equal to   ), that means that the angle between the two vectors is equal to  (=) 90 degrees
    If dot_prod < 0  (less than  ), that means that the angle between the two vectors is greater than (>) 90 degr.
    The dot product always takes the center of the grid in consideration to calculate the angles.

    If we use unit vectors to calculate the dot product, we get an interesting result!
        If they both face the same direction (parallel 0°), the scalar is 1
        If they face towards opposite direction (parallel 180°), the scalar is -1

        This is cool, because it is exactly the angle of a cosine function.

        !!! angle_in_radians = np.cos ( np.dot(a,b)) !!!

    PLANE
    Imagine that perpendicular to the vector (and through the origin), there passes a plane. 
    Unit normal vectors are unit vectors that describe the direction of the surface (de normaalhoek, want die staat loodrecht op het vlak). 

    Distance to that plane: 
        The dot product between a unit vector and any point in space( dot product between vector and position), returns the distance from the point to the plane
        In Godot language, there is the dotproduct function that returns " distance = np.dot(normal_vector, point_in_space)
            So I'm guessing a point and a vector are different

        Let D = scalar, which is the distance from the origin. Don't confuse eigenvalue scalar with just a floating point scalar.
        N (normal) has its own starting point at a plane. That plane is different from either x-axis or y-axis. D is that distance from the origin of the normal to the origin of the cartesian system.
"""

import numpy as np
import json


class Nucleoside:

    def __init__(self,jsonfile):
        self.splitted = jsonfile.split('.')[0]

        with open(jsonfile, 'r') as jsonf:
            self.jason = json.load(jsonf)

    #
    # Create functions to calculate for dihedrals
    def get_array(self):

        return np.asarray(json.loads(self.jason['pdb_properties']['Coordinates']), dtype=float)



class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """
    pass

