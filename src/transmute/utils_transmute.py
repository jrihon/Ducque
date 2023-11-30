import sys
from os.path import isfile, basename, dirname
from os import getcwd
from random import randint
import numpy as np

import systemsDucque as SD
from ducquelib.library import TABLE_CHEMISTRY, TABLE_LINKER, TABLE_NUCLEOBASE

class TransmuteToJson:
    """ This class exists to convert any *.pdb type format to an appriopriate *.json type format.

        What is done here is that the pdb file prompted is parsed for all the necessary data to be able to build it back up again.

        --pdb `*.pdb read from`
        --chemistry `residue name`
        --conformation `pyranose or furanose or phi-psi conformation`
        --moietyType `nucleoside / linker`
        --bondangles `see Ducque manual`
        --dihedrals `see Ducque manual`

        The other flags are easily parsed from the input file.  """


    def __init__(self, pdbfile):
        """ Initialise the object and create object properties"""
#        DUCQUEHOME = SD.return_DUCQUEHOME()
#        JSONDIR = DUCQUEHOME + "/json"

        self.fileName = pdbfile                     # path to the pdb file we input
        self.array = np.array([])
        self.atomName = list()
        self.residueName = str
        self.elementSymbol = list()


    def pdb_for_attributes(self):
        """ Reads the name of the file and converts the entire file into a workable dataframe.
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
            Charge:           line 79 - 80          """

        # Start new lists to append it all
        atomName, resName, xCoords, yCoords, zCoords, elementSymbol = ([] for _ in range(6))

        # Check if file is in cwd or in the pdb directory
        pdbfname = self.fileName
        if not isfile(pdbfname) :
            SD.print_filenotfound(pdbfname)
            SD.exit_Ducque()

        # Read the file and fill out the dataframe
        with open(pdbfname) as pdbfile:
            for line in pdbfile:
                if line[:4] == 'ATOM' or line[:6] == 'HETATM':

                    _name = line[12:16]
                    atomName.append(_name)

                    _res = line[17:20]
                    resName.append(_res)

                    _xcoord = line[30:38]
                    xCoords.append(_xcoord)

                    _ycoord = line[38:46]
                    yCoords.append(_ycoord)

                    _zcoord = line[46:54]
                    zCoords.append(_zcoord)

                    elemSym = line[76:78]
                    elementSymbol.append(elemSym)

            self.atomName = atomName
            self.residueName = resName[0].strip()
            self.elementSymbol = elementSymbol

            # Add the array as an attribute
            self.array = np.array([xCoords, yCoords, zCoords], dtype=float).T

    def get_array(self) -> list:
        """ Create a list of the array of coordinates respectively."""
        return np.ndarray.tolist(self.array)


    def get_shape_array(self) -> list:
        """ Get the shape of the array, but in a list instead of a tuple"""
        return [self.array.shape[0], self.array.shape[1]]


    def get_atoms(self) -> list:
        """ Get atom names to a list, also strip any remaining whitespace in all the strings in the list."""
        return list(map(lambda x: x.strip(), self.atomName))


    def get_element_symbol(self):
        """ Get element symbol to a list, also strip any remaining whitespace in all the strings in the list."""
        return list(map(lambda x: x.strip(), self.elementSymbol))


    def get_chemistry(self) -> str:
        """ Get the chemistry of the nucleoside, also strip any remaining whitespace in all the strings in the list. """
        return self.residueName


    def get_full_name(self, chemistry : str, moietyType) -> str:
        """ Get the full name of the nucleic acid chemistry or linker moietyType we want to convert to a json """

        if moietyType.upper() == "NUCLEOSIDE":
            try : 
                chemistry = TABLE_CHEMISTRY[chemistry.upper()]
            except :
                SD.print_invalid_key(chemistry, "TABLE_CHEMISTRY")
                sys.exit(1)
            
            return chemistry

        if moietyType.upper() == "LINKER":
            try : 
                 link = TABLE_LINKER[chemistry.upper()]
            except :
                SD.print_invalid_key(chemistry.upper(), "TABLE_LINKER")
                sys.exit(1)
            
            return link


    def get_nucleobase(self, nucleobase: str) -> str:
        """ Get the base that corresponds with this nucleic acid. Take the last character of the string Residue Name and
            look for it in the dictionary """
        try : 
            base = TABLE_NUCLEOBASE[nucleobase.upper()]
        except :
            SD.print_invalid_key(nucleobase, "TABLE_NUCLEOBASE")
            sys.exit(1)
        
        return base


    def get_angles(self, moietyType : str, angles_list : list) -> dict:
        """
        List the bond angles differently whether it belongs to a nucleoside or a linker moietyType 
        This function is used by both the bondangle parser and the dihedral parser.

        """
        json_dict = {}

        if moietyType == "nucleoside":
            # Strip the list of (for now) string values of their comma 
            angles_list = list(map(lambda x: x.strip(","), angles_list))
            if len(angles_list) == 7:
                angles_of_interest = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]
            elif len(angles_list) == 8:
                angles_of_interest = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "nu", "chi"]
            elif len(angles_list) == 6:
                angles_of_interest = ["alpha", "beta", "gamma", "delta", "epsilon", "chi"]
            else:
                print("Amount of dihedrals prompted is not aligned with the standard amount of dihedrals.")
                SD.exit_Ducque()


            # Check if all values are float
            for i, val in enumerate(angles_list):
                if not isinstance(float(val), float):
                    SD.print_conversion_err(angles_of_interest[i], val)
                    SD.exit_Ducque()

            # Check if size of the prompted dihedral values is the same as the amount of required dihedrals
            assert len(angles_of_interest) == len(angles_list),  "Check your input for missing dihedral values or missplaced commas.\nNote: the decimal values should be denoted by a point and not a comma."

            # Append the values to their respective bond angles
            for ang in range(len(angles_list)):
                json_dict[angles_of_interest[ang]] = float(angles_list[ang])

            return json_dict 


        elif moietyType == "linker":
            angles_list = list(map(lambda x: x.strip(","), angles_list))

            json_dict = {}

            # Check if all values are float
            for i, val in enumerate(angles_list):
                if not isinstance(float(val), float):
                    SD.print_conversion_err(angles_list[i], val)
                    SD.exit_Ducque()

            # Bond Angles
            if len(angles_list) == 1:
                json_dict["angle_1"] = float(angles_list[0])

            # Dihedrals
            elif len(angles_list) == 2:
                json_dict["dihedral_1"] = float(angles_list[0])
                json_dict["dihedral_2"] = float(angles_list[1])

            return json_dict 

        else : 
            SD.print_invalid_argument(moietyType, "--moiety")


    def get_output_name(self, chemistry : str, moietyType : str, conformation : str, nucleobase: str) -> str:
        """ Create the name of the file based on the chemistry of the nucleic acid chemistry and its corresponding base 
            This function creates the name of the json file

            If moietyType == linker, then the nucleobase == "" (empty). 
            The conformation is optional and is prompted as either R or S stereochemistry.
        """


        if moietyType == "nucleoside":
            name_of_chemistry = chemistry.lower()
            name_of_base = self.get_nucleobase(nucleobase).lower()

            conformation = conformation.lower()
            return name_of_chemistry + "_" + name_of_base + "_" + conformation

        elif moietyType == "linker":
            name_of_chemistry = chemistry.lower()
            name_of_linker = TABLE_LINKER[chemistry].lower()
            
            if conformation.upper() == "R" or conformation.upper() == "S":
                return  name_of_chemistry + "_" + name_of_linker + "_" + conformation.lower()
            elif conformation.upper() != "NONE":
                SD.print_invalid_argument(conformation, "--conformation")
            else:
                return name_of_chemistry + "_" + name_of_linker

        else :
            sys.exit("The molecule is not annotated with either `nucleoside` or `linker`. Please revise the inputs")


class TransmuteToPdb:
    """ This class is used to convert the *.xyz files from ORCA to *.pdb files. Later on, these *.pdb files are prompted into Ducque to convert to *.json files.

        --xyz `*.xyz`
        --residue `Residue Name`
        --atomname_list `see Ducque manual`               """


    def __init__(self, xyzfile):
        """ Initialise the object and create object attributes """
        self.pathName = xyzfile
        self.rootName = basename(xyzfile)[:-4] # cut off last `.xyz` part
        self.pdbName = self.rootName + ".pdb"
        self.fullPathTo = getcwd() + "/" + dirname(xyzfile)
        self.x = list()
        self.y = list()
        self.z = list()
        self.elements = list()
        self.atomNameList = list()


    def parse_xyz_and_elementsymbol(self):
        """ Reads the xyz datafile and returns the coordinates and the element symbol"""
        # read only the lines with x-y-z coordinates
        if not isfile(self.pathName) :
            SD.print_filenotfound(self.pathName)
            SD.exit_Ducque()

        with open(self.pathName, "r") as XYZ :
            _fileList = [line.strip() for line in XYZ.readlines()[2:]]

        # extract the coordinates. 
        xCoords, yCoords, zCoords, elements = [], [], [], []
        for line in _fileList:
            x, y, z = '{:.3f}'.format(float(line.split()[1])), '{:.3f}'.format(float(line.split()[2])), '{:.3f}'.format(float(line.split()[3]))
            ele = line.split()[0].strip()

            xCoords.append(x)
            yCoords.append(y)
            zCoords.append(z)
            elements.append(ele)

        self.x = xCoords
        self.y = yCoords
        self.z = zCoords
        self.elements = elements


    def return_processed_atomname_list(self, atomNameList : list):
        """ Returns a list of atom names in a well formatted list """
        atomNameList = list(map(lambda x : x.strip(","), atomNameList))

        for atom in atomNameList:
            if len(atom) > 4:
                print(f"The following atom has too many characters in the string {atom}. Maximum amount allowed is 4.\n")
                SD.exit_Ducque()

        self.atomNameList = atomNameList


    def arraysize_vs_atomname_list_compatibility(self) -> bool:
        """ Checks to see if the size of the parsed cartesian coordinates (in length) matches the size of the inputted atomname_list. They should match! """

        if len(self.elements) == len(self.atomNameList):
            return True
        else :
            print(f"Atomname_list length : {len(self.atomNameList)}.\nCoordinate array size : ({len(self.elements)}, , 3 ).")
            return False


    def elementsymbol_vs_atomname_list_compatibility(self) -> bool:
        """ Checks if the same type of atoms are prompted in atomname, when comparing them to the element symbol list """

        list_of_two_character_valid_atoms = ["Cl", "Br", "Si"]

        for i_atom in range(len(self.atomNameList)):

            # Check if the atomname you are testing is longer than two characters. If it is, check if it belongs to a valid element in the PSE
            if self.elements[i_atom] in list_of_two_character_valid_atoms:
                atom_type = self.atomNameList[i_atom][:2]

                # if it does not belong to a valid element, either add to the list up here or change the element to something that exists
                if atom_type != self.elements[i_atom]:
                    print(f"Check for capitalization of the prompted atom : {atom_type}.\n"
                            "If that did not prompt the error, check for its validity as an atom.\n"
                            "Here is the valid atom list : {list_of_two_character_valid_atoms}.\n"
                            "Feel free to adjust this to your needs in the file : Ducque/src/transmute/utils_transmute.py ; elementsymbol_vs_atomname_list_compatibility() function.\n")
                    return False

            # First parse out only the first character of the atomNameList, since this the first letter denotes the element of that atom
            atom_type = self.atomNameList[i_atom][0]

            if atom_type != self.elements[i_atom]:
                print(f"Atoms that do not match : Atom prompted - {atom_type}. Element parsed : {self.elements[i_atom]}. Position {str(i_atom + 1)}\n")
                return False

        return True


    def write_to_pdb_format_file(self, residue) :

        randomised_integer_for_sequence_number = randint(1,100)

        atomNumbers = np.linspace(1, len(self.elements), len(self.elements), dtype=int)

        # Write out the `*.pdb` file
        write_to_file = self.fullPathTo + "/" + self.pdbName
        with open(write_to_file, "w") as pdb:

            for idx in range(len(self.elements)):
                line = ["ATOM",
                        atomNumbers[idx],
                        self.atomNameList[idx],
                        residue,
                        "A",
                        randomised_integer_for_sequence_number,
                        self.x[idx], self.y[idx],self.z[idx],
                        "1.00", "0.00",
                        self.elements[idx] ]

                pdb.write("%-4s  %5d %-4s %3s %s%4d    %8s%8s%8s%6s%6s          %2s\n" % tuple(line))

        SD.print_writing(write_to_file)
