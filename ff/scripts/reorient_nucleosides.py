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


How this script works :
    Make sure that you have all the pdb files of your respective nucleosides, like so

    CWD
    |-dA      ---    conf1.pdb conf2.pdb ...
    |-dC      ---    conf1.pdb conf2.pdb ...
    |-dG      ---    conf1.pdb conf2.pdb ...
    |-dT      ---    conf1.pdb conf2.pdb ...
    |-linker  ---    linker.pdb

    The script will go through all the subdirectories in the current working directory (so make sure there arent more of them than just the nucleosides),
        parse all the pdb files in the respective subdirectories and then orient them accordingly.


ARGUMENTS :
    The scripts takes in two argument.

    MOLECULE_ARG : 'nucleoside' or 'linker' when you want to reorient the linker
        This will let the script know whether we want to reorient a group of nucleosides or just a single linker moiety

"""

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
            print(f"Could not find {pdbfname} in the directory.\n")

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




def arrange_pdbfiles_per_nucleoside(listPdb : list) -> str:
    """ Parse the correct subdirectory to then add the string to the name of the file we are going to generate """
    # Be sure at all files are in the correct subdir
    subdir = [x.split("/")[0] for x in listPdb ][0]

    return subdir


def sort_by_integer(listPdb : list) -> list:
    """ Sort the filenames in the list by their given order
        example : [MinOpt02.pdb, MinOpt03.pdb, MinOpt01.pdb] """

    # Assumes that you have 2 integers denoting the number of the pdb files ; MinOpt02.pdb (instead of MinOpt2.pdb)
    int_list = [ int(x.split(".")[0][-1]) for x in listPdb]

    sorted_int_list = sorted(int_list)

    # Retrieve the indexes from the values that match with the new list, so we can create a new list with the correct order
    idx_list = np.empty(shape=len(int_list), dtype=int)
    for i, _int in enumerate(int_list):
        idx_list[i] = sorted_int_list.index(_int)

    # final list
    final_list = list(range(len(int_list)))
    for i in range(len(int_list)):
        final_list[idx_list[i]] = listPdb[i]

    return final_list


def return_all_pdb_files_in_subdirectories(cmd_arg : str) -> list:
    """ Go through all the subdirectories in the os.get_cwd() and retrieve all the pdb files in there.
        The os.get_cwd() should be so that all the required pdb files should be in their respective subdirectories like :
            cwd/ -> Adenine/, Guanine/, Cytosine/, Thymidine/   """

    # current directory
    cwd = os.listdir(os.getcwd())

    allPdbFiles = list()

    # Searching for the linker directory
    if cmd_arg == "linker":
        for _dir in cwd :
            if os.path.isdir(_dir) :

                # get all pdb files
                if "LINKER" in _dir.upper() :
                    listContent = [_dir + "/" + x for x in os.listdir(_dir) if x.endswith(".pdb")]
                else:
                    continue

                # if no pdb files in this directory
                if len(listContent) == 0:
                    sys.exit(f"The {_dir} directory contains no pdb files\n")

                for pdbfile in listContent:
                    if "_model.pdb" in pdbfile:
                        print(f"The {_dir} directory already contains a *_models.pdb\n")
                        sys.exit(0)
                allPdbFiles.append(listContent[0])
                return allPdbFiles


    # If there are any subdirectories, get the listed content and only get the pdb files from it
    for _dir in cwd :
        if os.path.isdir(_dir) :
            if "LINKER" in _dir.upper() :
                continue

            pdbFOUND = False

            # get all pdb files
            listContent = [_dir + "/" + x for x in os.listdir(_dir) if x.endswith(".pdb")]

            # if no pdb files in this directory
            if len(listContent) == 0:
                print(f"The {_dir} directory contains no pdb files\n")
                continue

            # Check if the reorient_nucleosides.py has already been used
            for pdbfile in listContent:
                if "_model.pdb" in pdbfile:
                    print(f"The {_dir} directory already contains a *_models.pdb\n")
                    pdbFOUND = True
                    continue

            if not pdbFOUND:
                listContent = sort_by_integer(listContent)

                # Append instead of extend, because it makes it easier to parse the correct files from the same directory
                allPdbFiles.append(listContent)

    return allPdbFiles


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




#-------------------------- MAIN --------------------------
def main():

    # Check validity of first prompted argument
    try:
        MOLECULE_ARG = sys.argv[1]
        if not MOLECULE_ARG in ["nucleoside", "linker"]:
            sys.exit(f"This argument does not take in : {MOLECULE_ARG}\n"
                    "Possibility : \n"
                    "   nucleoside\n"
                    "   linker\n")
    except IndexError:
        sys.exit("No argument for the type of molecule has been prompted for MOLECULE_ARG\n")


    ### Start program

    if MOLECULE_ARG == "nucleoside":
        # NUCLEOSIDE
        atoms = ["C1'", "C5'", "N3'"]
        atomOne = atoms[0]
        atomTwo = atoms[1]
        atomThree = atoms[2]

        # Create a reference nucleoside where we superimpose all the other nucleosides on
        allPdbFiles = return_all_pdb_files_in_subdirectories(MOLECULE_ARG)
        REFERENCE_NUCLEOSIDE = allPdbFiles[0][0]
        pdb_instance = pdbMolecule(REFERENCE_NUCLEOSIDE)



    if MOLECULE_ARG == "linker":
        # LINKER
        atoms = ["P", "OP1", "OP2"]
        atomOne = atoms[0]
        atomTwo = atoms[1]
        atomThree = atoms[2]

        allPdbFiles = return_all_pdb_files_in_subdirectories(MOLECULE_ARG)
        REFERENCE_NUCLEOSIDE = allPdbFiles[0]
        pdb_instance = pdbMolecule(REFERENCE_NUCLEOSIDE)



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

    # Let's remember the C1' - > N3' vector so we can rotate onto it later
    vec2 = return_normalized(nuc_array[array_of_atoms[2]])
    vec1 = return_normalized(nuc_array[array_of_atoms[1]])
    directionAxis1 = get_direction_of_rotation(vec1, vec2)




    # LINKER
    if MOLECULE_ARG == "linker":
        _dir = REFERENCE_NUCLEOSIDE.split("/")[0]
        with open(_dir + "/LINKER_model.pdb" , "w") as pdb:
            pdb_instance.Write_Out_Pdb(nuc_array, pdb)

            pdb.write("END")
            sys.exit(0)



    # NUCLEOSIDE
    if MOLECULE_ARG == "nucleoside":
        # Iterate over all the files
        for i, nucleobase in enumerate(allPdbFiles):

            # Create a new model file
            if MOLECULE_ARG == "nucleoside":
                subdir = arrange_pdbfiles_per_nucleoside(nucleobase)

            # Write out the models
            for j, conformation in enumerate(nucleobase):

                with open("./" + subdir + "/" + subdir + "_model" + str(j + 1) + ".pdb" , "w") as pdb:
#                    pdb.write("MODEL " + str(j + 1) + "\n")
                    pdb_instance = pdbMolecule(conformation)
                    pdb_instance.pdb_to_dataframe()

                    array_of_atoms = pdb_instance.retrieve_atom_index_MULTIPLE([atomOne, atomTwo, atomThree])

                    # Get nucleoside array to origin, also instantiate a duplicate array to work with
                    nuc_array = pdb_instance.array - pdb_instance.array[array_of_atoms[0]]

                    # Turn the first time onto the x-axis [1,0,0]
                    vAxis1 = np.array([1,0,0])
                    vAxis2 = return_normalized(nuc_array[array_of_atoms[1]])
                    quaternion1 = get_quaternion(vAxis1, vAxis2)
                    nuc_array = quaternion1.apply(nuc_array)

                    # Second rotation on the y-axis [0,1,0]
                    vec4 = return_normalized(nuc_array[array_of_atoms[2]])
                    vec3 = return_normalized(nuc_array[array_of_atoms[1]])
                    directionAxis2 = get_direction_of_rotation(vec3, vec4)
                    quaternion2 = get_quaternion(directionAxis1, directionAxis2)
                    nuc_array = quaternion2.apply(nuc_array)

                    # Write out the pdb as a test
                    pdb_instance.Write_Out_Pdb(nuc_array, pdb)

                    pdb.write("END")


if __name__ == "__main__":
    main()


