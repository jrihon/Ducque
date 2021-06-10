import numpy as np
import pandas as pd
import json, sys

import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2


"""
The functions that run the matrix rotations
"""

# CLASSES
class Nucleoside:

    def __init__(self, jsonfile):

        with open(jsonfile, "r") as jsonf:
            self.jason = json.load(jsonf)

        self.splitted = jsonfile.split(".")[0]
        self.array =  np.asarray(json.loads(self.jason["pdb_properties"]["Coordinates"]), dtype=float)
        self.atom_list = json.loads(self.jason["pdb_properties"]["Atoms"])


    def get_alpha(self) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["alpha"])

    def get_beta(self) -> float:
        # because beta is still inside a the backbone dictionary, we need to load it again 
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["beta"])

    def get_zeta(self) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["zeta"])

    def get_epsilon(self) -> float:
        return float(json.loads(self.jason["Dihedrals"]["Backbone"])["epsilon"])


class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_COP(self) -> float:
        # since the angles are inside the angles dictionary, we don"t need to load string again
        return float(self.jason["Angles"]["C5_O5_P"]) * (np.pi / 180)

    def get_P(self) -> float:
        return self.array[0]

    def get_OPO2(self) -> float:
        return float(self.jason["Angles"]["O5_P_OP2"]) * (np.pi / 180)

    def get_OPO1(self) -> float:
        return float(self.jason["Angles"]["O5_P_OP1"]) * (np.pi / 180)

    def get_O3PO5(self) -> float:
        return float(self.jason["Angles"]["O5_P_O3"]) * (np.pi / 180)

    def get_OP2_dihedral(self) -> float:
        return float(self.jason["Dihedrals"]["dihedral_oxygen_OP2"])

    def get_OP1_dihedral(self) -> float:
        return float(self.jason["Dihedrals"]["dihedral_oxygen_OP1"])


## FUNCTIONS THAT ARE MEANT TO BYPASS THE ITERATIVE CODING AND RESULT IN ONLY THE NECESSARY RESULTS IN THE LABYRINTH.PY : ARCHITECTURE()
def generate_vector_of_interest(angle : float, dihedral : float, atom_array : np.array) -> np.array:
    """ This function generates a single vector for which there exists only one angle and dihedral.
    The atom_array contains the three first atoms in the sequence that make up the dihedral.
    Example = if the sequence of a dihedral is C4' - C5' - O5' - P (beta backbone), then the atom_array is [O5', C5', C4']
    """
    #print(1)
    # Get the vector to rotate the cone vector onto. This is done on the middle two atoms of the sequence
    # Example : [O5'] minus (-) [C5'] results in a vector that goes [C5' -> O5']
    p0 = LFT1.return_normalized(atom_array[0] - atom_array[1])

    # Generate the cone vector
    cone_vector = LFT1.generate_cone_vector(angle)

    # Get the quaternion that corresponds to the desired rotation<
    # Note, if 'vector_to_from" is empty, it defaults to the Z-axis, around which the cone is generated np.array([0,0,1])
    quaternion = LFT1.get_quaternion(p0 * -1.0)

    # Rotate the cone vector
    rotated_cone_vector = LFT1.rotate_with_quaternion(quaternion, cone_vector)

    # Calculate the possible dihedrals we have with the first three atoms and the cone vector. Remember: the cone vector represents the position of the fourth atom
    range_of_dihedrals = LFT1.dihedral_array(atom_array, rotated_cone_vector)

    # Now extrapolate in between the range_of_dihedrals to get the correct theta angle to generate the single vector of interest.
    # The theta angle corresponds to the angle (spherical coordinates r, theta, phi) at which we need to generate our vector at.
    theta_interpolate = LFT1.get_interpolated_dihedral(range_of_dihedrals, dihedral)

    # Generate the correct vector of interest
    single_vector = LFT1.generate_and_rotate_single_vector_QUAT(theta_interpolate, angle, quaternion)

    return single_vector


def create_PDB_from_matrix(matrix : np.ndarray, list_of_sequence : list) -> None:
    print("Writing to pdb ...")

    # Parse From json
    #   -AtomName : CHECK
    #   -Sequence should be deduced from the shape of the nucleoside and the linker : CHECK
    #   -Residue name, but can't be arsedi : CHECK
    #   - Elementsymbol parsed the sane way that AtomName will be parsed : CHECK
    df = pd.DataFrame()

    df['RecName'] = ['ATOM' for x in range(matrix.shape[0])]
    df['AtomNum'] = np.arange(start=1, stop=matrix.shape[0] + 1)
    df['AtomName'] = LFT2.pdb_AtomNames_or_ElementSymbol(list_of_sequence, "Atoms")
    df['AltLoc'] = ' '
    df['ResName'] = LFT2.pdb_Residuename(list_of_sequence)
    df['Chain'] = 'A'
    df['Sequence'] = LFT2.pdb_Sequence(list_of_sequence)
    df['X_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,0]))
    df['Y_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,1]))
    df['Z_coord'] = list(map(lambda x: '{:.3f}'.format(x), matrix[:,2]))
    df['Occupancy'] = '1.00'
    df['Temp'] = '0.00'
    df['SegmentID'] = str('   ')
    df['ElementSymbol'] = LFT2.pdb_AtomNames_or_ElementSymbol(list_of_sequence, "Symbol")

    with open('testing_daedalus.pdb' ,'w') as pdb:
        for index, row in df.iterrows():
            split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
            pdb.write('%-6s%5s%5s%s%s%3s%5s  %8s%8s%8s%6s%6s%4s      %2s\n' % tuple(split_line))
        pdb.write('END')

