import os
import sys
import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform.rotation import Rotation as R
from pandas import DataFrame

# LIBRARIES : NUMPY, SCIPY, PANDAS

# The first this strategy, the ﬁrst selected atom is translated to the origin of axes, the ﬁrst two atoms deﬁne the (0, X) axis while the third one is used to
# deﬁne the (0, X, Y) plane. The (0, Z) axis is automatically set as the cross-product between the (0, X) and (0, Y) axes. This approach can be used for every optimized molecular
# geometry, and is the basis for multiple orientation charge ﬁtting.

#---- Pseudo-code
#   1. Define which atoms you want to rotate on.
#       - The C1' will be at the center of the cartesian system; C1' location([0,0,0])
#       - The C6' will be the axis([1,0,0]) on which to impose the first rotation (C1' -> C6')
#       - The N3' will be the next axis([0,0,1]) as the direction axis, so the three selected atoms will lie in the 0-X-Y plane.
#       - This plane will be the plane upon which all other nucleosides are rotated on, reference from the first nucleoside that is being parsed.
#   


"""
This scripts requires three libraries in python : NumPy, SciPy, Pandas

IMPORTANT : in the main() function there is a variable called atoms.
    -- 'atoms' is a list variable, where the names of the atoms are contained in.
        changes this to your needs. Remember that the three atoms are in a plane.
        If possible, make it consistent with how DNA/RNA are superposed onto the (0, X, Y).



ARGUMENTS :
    The scripts takes in one argument, the pdb file

"""

#-------------------------- MAIN --------------------------
def main():

    # Check validity of first prompted argument
    try:
        MOLECULE_ARG = sys.argv[1]
    except IndexError:
        print("No arguments prompted")

    if not os.path.isfile(MOLECULE_ARG):
        sys.exit(f"FileNotFoundError : {MOLECULE_ARG} is not present. ")


    ### Start program
    atoms = ["C1'", "C2'", "C4'"]
    atomOne = atoms[0]
    atomTwo = atoms[1]
    atomThree = atoms[2]

    pdb_instance = pdbMolecule(MOLECULE_ARG)


    pdb_instance.pdb_to_dataframe()
    array_of_atoms = pdb_instance.retrieve_atom_index_MULTIPLE([atomOne, atomTwo, atomThree])

    # Get nucleoside array to origin, also instantiate a duplicate array to work with
    nuc_array = pdb_instance.array - pdb_instance.array[array_of_atoms[0]]

    # Turn the first time onto the x-axis [1,0,0]
    Xaxis = np.array([1,0,0])
    vec0 = return_normalized(nuc_array[array_of_atoms[1]])
    quatX = get_quaternion(Xaxis, vec0)
    nuc_array = quatX.apply(nuc_array)

    # Turn the second time onto the X-Z plane, by rotating the direction axis onto the y-axis [0,0,1]
    vec2 = return_normalized(nuc_array[array_of_atoms[2]])
    vec1 = return_normalized(nuc_array[array_of_atoms[1]])
    directionAxis1 = get_direction_of_rotation(vec1, vec2)
    Yaxis = np.array([0,0,1])
    quatY = get_quaternion(Yaxis, directionAxis1)
    nuc_array = quatY.apply(nuc_array)


    with open("./reoriented_" + pdb_instance.filename, "w") as pdb:
        pdb_instance.Write_Out_Pdb(nuc_array, pdb)

 
class pdbMolecule:

    def __init__(self, pdbfile):
        """ Initialise the object and create object properties"""
        self.splitted = pdbfile.split('.')[0]
        self.pdb_dataframe = DataFrame()
        self.array = np.array([])
        self.filename = pdbfile

    def pdb_to_dataframe(self):
        """
        Reads the name of the file and converts the entire file into a workable dataframe.
        By convention of the columns, below, the dataframe is set up like this.
        __________________________________________________________________________________________________________
            # Clarification of the column headers
        Record name :     line 1 - 6
        Atom number :     line 7 - 11
        Atom name :       line 13 - 16
        alternat. loc.:   line 17      (Avaline and Bvaline, delete Bvaline and then clear the column)
        Residue name:     line 18 - 20
        Chain :           line 22
        Sequence number:  line 23 - 26
        X_coord:          line 31 - 38
        Y_coord:          line 39 - 46
        Z_coord:          line 47 - 54
        Occupancy:        line 55 - 60 (no idea what this is; keep)
        Temp. factor:     line 61 - 66 (no idea what this is either; keep)
        Segment id:       line 73 - 76 (no idea what this is, keep it and clear column)
        Element symbol:   line 77 - 78
        Charge:           line 79 - 80
        """
        # Start new lists to append it all
        AtomName, ResName, X_coords, Y_coords, Z_coords, ElementSymbol = ([] for i in range(6))

        # Check if file is in cwd or in the pdb directory
        pdb_fname = self.filename
        try:
            os.path.isfile(pdb_fname)
        except FileNotFoundError:
            print(f"Could not find {pdb_fname} in the directory.\n")

        # Read the file and fill out the dataframe
        with open(pdb_fname) as pdbfile:
            for line in pdbfile:
                if line[:4] == 'ATOM' or line[:6] == 'HETATM':

                    name = line[12:16]
                    AtomName.append(name)

                    resname = line[17:20]
                    ResName.append(resname)

                    _Xcoord = line[30:38]
                    X_coords.append(_Xcoord)

                    _Ycoord = line[38:46]
                    Y_coords.append(_Ycoord)

                    _Zcoord = line[46:54]
                    Z_coords.append(_Zcoord)

                    ElemSym = line[76:78]
                    ElementSymbol.append(ElemSym)

            self.pdb_dataframe['AtomName'] = AtomName
            self.pdb_dataframe['ResName'] = ResName
            self.pdb_dataframe['X_Coord'] = X_coords
            self.pdb_dataframe['Y_Coord'] = Y_coords
            self.pdb_dataframe['Z_Coord'] = Z_coords
            self.pdb_dataframe['ElementSymbol'] = ElementSymbol

            # Add the array as an attribute
            self.array = np.array([X_coords, Y_coords, Z_coords], dtype=float).T

            # Add the atom name list as an attribute
            self.atom_list = list(map(lambda x : x.strip(), AtomName))

            # Add elements as an attribute
            self.elements = list(map(lambda x : x.strip(), ElementSymbol))

            # Add residue name
            self.resname = ResName[0].strip()


    def retrieve_atom_index_MULTIPLE(self, atoms : list, index_counter : int = 0) -> np.array :
        """ Retrieves the index in the self.array of the atom of interest
            This integer will be used to retrieve the vector of the atom of interest """
        array_of_indexes = np.zeros(len(atoms), dtype=int)

        for i in range(len(atoms)):
            array_of_indexes[i] = self.atom_list.index(atoms[i]) + index_counter

        return array_of_indexes


    def Write_Out_Pdb(self, atom_array : np.ndarray, pdb):
        """ Write out the data for the pdb filename
            @param : pdb . This is a file that has already opened
            https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html"""

        # LEADING STRAND
        df_nucleoside = DataFrame()

        df_nucleoside["RecName"] = ["ATOM" for x in range(atom_array.shape[0])]
        df_nucleoside["AtomNum"] = np.arange(start=1, stop=atom_array.shape[0] + 1)
        df_nucleoside["AtomName"] = self.atom_list
        df_nucleoside["AltLoc"] = " "
        df_nucleoside["ResName"] = self.resname
        df_nucleoside["Chain"] = "J"
        df_nucleoside["Sequence"] = 1
        df_nucleoside["X_coord"] = list(map(lambda x: "{:.3f}".format(x), atom_array[:,0]))
        df_nucleoside["Y_coord"] = list(map(lambda x: "{:.3f}".format(x), atom_array[:,1]))
        df_nucleoside["Z_coord"] = list(map(lambda x: "{:.3f}".format(x), atom_array[:,2]))
        df_nucleoside["Occupancy"] = "1.00"
        df_nucleoside["Temp"] = "0.00"
        df_nucleoside["SegmentID"] = str("   ")
        df_nucleoside["ElementSymbol"] = self.elements

        # Write out the pdb file
        for index, row in df_nucleoside.iterrows():
            split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
            pdb.write("%-6s%5s%5s%s%3s%2s%5d  %8s%8s%9s%6s%7s%4s     %2s\n" % tuple(split_line))



# MATHEMATICS
def return_normalized(vector : np.ndarray) -> np.ndarray:
    """ returns a normalized vector """
    return vector / norm(vector)


# MATHEMATICS
def get_direction_of_rotation(from_vector : np.ndarray, vector_to_rotate_onto : np.ndarray) -> np.ndarray:
    """ Cross product with the cone's axis to get the direction of the rotation axis
        Get the direction where vectorA rotates onto vectorB ; u = vectora X vectorb
        Let's normalize the direction """

    cross_product = np.cross(from_vector, vector_to_rotate_onto)
    return return_normalized(cross_product)


# MATHEMATICS
def get_angle_of_rotation(from_vector : np.ndarray, vector_to_rotate_onto : np.ndarray) -> float:
    """ The scalar product (dot product) to get the cosine angle.
        Here we do the arccos, to get the angle immediately. """
    return np.arccos(np.dot(from_vector, vector_to_rotate_onto))


# MATHEMATICS
def get_quaternion(vector_to_rotate_onto : np.ndarray, vector_to_rotate_from : np.ndarray ) -> np.array:
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




if __name__ == "__main__":
    main()
