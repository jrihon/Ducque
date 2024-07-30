from io import TextIOWrapper
import os
import sys
import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation

# DEPENDECIES : NUMPY, SCIPY

# The first this strategy, the ﬁrst selected atom is translated to the origin of axes, the first two atoms define the X-axis while the third one is used to
# define the (X, Y) plane. The Z-axis is the normal to the cross-product between the X- and Y-axis. This approach can be used for every optimized molecular
# geometry, and is the basis for multiple orientation charge fitting.

#---- Pseudo-code
#   1. Define which atoms you want to rotate on.
#       - Suppose `atoms = ["C1'", "C6'", "N3'"]`
#       - The C1' will be at the center of the cartesian system; C1' location([0,0,0])
#       - The C6' will be the axis([1,0,0]) on which to impose the first rotation (vector :: C1' -> C6')
#       - The N3' will be the next axis([0,0,1]) as the direction axis, so the three selected atoms will lie in the ORIGIN-X-Y plane.
#       - This plane will be the plane upon which all other nucleosides are rotated on, reference from the first nucleoside that is being parsed.
#   


"""
This scripts requires the following dependecies in Python3 : NumPy, SciPy

IMPORTANT : in the main() function there is a variable called atoms.
    -- 'atoms' is a list variable, where the names of the atoms are contained in.
        changes this to your needs. Remember that three atoms form a plane.



ARGUMENTS :
    The scripts takes in one argument, the molecule's pdb file in need of reorienting

"""

#-------------------------- MAIN --------------------------
def main():

    # Check validity of first prompted argument
    try:
        MOLECULE_ARG = sys.argv[1]
    except IndexError:
        print("No arguments prompted. Prompt molecule pdb file to reorient.")
        exit(1)

    if not os.path.isfile(MOLECULE_ARG):
        sys.exit(f"FileNotFoundError : {MOLECULE_ARG} is not present. ")


    #------------------------------
    #------------------------------
    #------------------------------
    # ADD ATOM NAMES OF PLANE MOLECULES HERE
    atoms = ["C1'", "C2'", "C4'"]

    atomOne = atoms[0]
    atomTwo = atoms[1]
    atomThree = atoms[2]
    #------------------------------
    #------------------------------
    #------------------------------




    # start program
    pdb_instance = pdbMolecule(MOLECULE_ARG)
    pdb_instance.parse_data()

    # get indices of atoms by name reference
    indices_of_atoms = pdb_instance.retrieve_atom_index_MULTIPLE([atomOne, atomTwo, atomThree])

    # move coordinate array to origin
    nuc_array = pdb_instance.coordinateArray - pdb_instance.coordinateArray[indices_of_atoms[0]]

    # Rotate the first time onto the x-axis [1,0,0]
    Xaxis = np.array([1,0,0])
    vec0 = return_normalized(nuc_array[indices_of_atoms[1]])
    quatX = get_quaternion(Xaxis, vec0)
    nuc_array = quatX.apply(nuc_array)

    # Rotate the second time onto the X-Y plane, by rotating the coordinate array's direction axis onto the z-axis [0,0,1]
    vec2 = return_normalized(nuc_array[indices_of_atoms[2]])
    vec1 = return_normalized(nuc_array[indices_of_atoms[1]])
    directionAxis1 = get_direction_of_rotation(vec1, vec2)
    Zaxis = np.array([0,0,1])
    quatY = get_quaternion(Zaxis, directionAxis1)

    # push reoriented coordinate array back into the pdbMolecule()'s attribute
    pdb_instance.coordinateArray = quatY.apply(nuc_array)


    # write to file
    with open("./reoriented_" + pdb_instance.filename, "w") as pdb:
        pdb_instance.write_out_pdb(pdb)


class PdbDataFrame:
    def __init__(self) -> None:
        pass
 
class pdbMolecule:

    def __init__(self, pdbfile):
        """ Initialise the object and create object properties"""
        self.splittedFilename = pdbfile.split('.')[0]
#        self.pdb_dataframe = PdbDataFrame()
#        self.array = np.array([])
        self.filename = pdbfile

    def parse_data(self):
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
        AtomName, ResName, X_coords, Y_coords, Z_coords, ElementSymbol = ([] for _ in range(6))

        # Check if file is in cwd or in the pdb directory
        pdb_fname = self.filename
        try:
            os.path.isfile(pdb_fname)
        except FileNotFoundError:
            print(f"Could not find {pdb_fname} in the directory.\n")
            exit(1)

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

            # Add the coordinate array
            self.coordinateArray = np.array([X_coords, Y_coords, Z_coords], dtype=float).T

            # Add the atom name list
            self.atom_list: list[str] = list(map(lambda x : x.strip(), AtomName))

            # Add element list
            self.elements = list(map(lambda x : x.strip(), ElementSymbol))

            # Add residue name
            self.resName = ResName[0].strip()


    def retrieve_atom_index_MULTIPLE(self, atoms : list) -> list :
        """ Retrieves the indices in the self.coordinateArray of the atom of interest
            These indices will be used to retrieve the vector of the atom of interest """

        return list(map(lambda atom: self.atom_list.index(atom), atoms))
#        indices_of_atoms = np.zeros(len(atoms), dtype=int)
#        for i in range(len(atoms)):
#            indices_of_atoms[i] = self.atom_list.index(atoms[i])
#
#        return indices_of_atoms


    def write_out_pdb(self, pdb: TextIOWrapper):
        """ Write out the data for the pdb filename
            https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html"""

        # precompute size
        sizeMolecule = len(self.atom_list)

        # Pdb Format File
        recordName = "ATOM"
        atomNumber = [x for x in range(1, sizeMolecule + 1)]
        # self.atom_list
        alternativeLocator = str(" ")
        # self.resname
        chainletter = "A"
        sequenceNumber = "1"
        #coordinates
        xCoord = list(map(lambda x: "{:.3f}".format(x), self.coordinateArray[:,0]))
        yCoord = list(map(lambda x: "{:.3f}".format(x), self.coordinateArray[:,1]))
        zCoord = list(map(lambda x: "{:.3f}".format(x), self.coordinateArray[:,2]))
        occupancy = "1.00"
        temperature = "0.00"
        segmentId = str("   ")
        # self.elements 

        # Write out the pdb file
        for idx in range(sizeMolecule):
            split_line = [ 
                          recordName, atomNumber[idx], self.atom_list[idx],
                          alternativeLocator, self.resName, chainletter, sequenceNumber,
                          xCoord[idx], yCoord[idx], zCoord[idx], 
                          occupancy, temperature, segmentId, self.elements[idx]
                          ]
            pdb.write("%-6s%5s%5s%s%3s%2s%5s  %8s%8s%9s%6s%7s%4s     %2s\n" % tuple(split_line))

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
def get_quaternion(vector_to_rotate_onto : np.ndarray, vector_to_rotate_from : np.ndarray ) -> Rotation:
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

    return Rotation.from_quat([qx, qy, qz, qw])




if __name__ == "__main__":
    main()
