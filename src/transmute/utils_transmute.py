import sys, os
import numpy as np
from typing import Union

import systemsDucque
import transmute.transmute_constants as TC

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
        DUCQUEHOME = systemsDucque.return_DUCQUEHOME()

        self.rootName = pdbfile.split('.')[0]
        self.fileName = DUCQUEHOME + "/" + self.rootName + ".pdb"
        self.array = np.array([])
        self.atomName = list()
        self.residueName = list()
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
        atomName, resName, xCoords, yCoords, zCoords, elementSymbol = ([] for i in range(6))

        # Check if file is in cwd or in the pdb directory
        pdbfname = self.fileName
        try:
            os.path.isfile(pdbfname)
        except FileNotFoundError:
            print(f"Could not find {pdbfname} in the directory.\n")

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
#
#            # Add the atom name list as an attribute
#            self.atomList = list(map(lambda x : x.strip(), AtomName))

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


    def get_full_name(self, identifier : str, moietyType):
        """ Get the full name of the nucleic acid chemistry or linker moietyType we want to convert to a json """

        if moietyType == "nucleoside":
            return TC.nucleoside_dict[identifier.upper()]

        if moietyType == "linker":
            return TC.linker_dict[identifier.upper()]


    def get_base(self) -> str:
        """ Get the base that corresponds with this nucleic acid. Take the last character of the string Residue Name and
            look for it in the dictionary """
        base = str(self.residueName)[-1]

        return TC.base_dict[base]


    def get_dihedrals(self, identifier : str, moietyType : str,  dihedrals_list : list) -> dict:
        """ List the dihedrals differently whether it belongs to a nucleoside or a linker moietyType """

        if moietyType == "nucleoside":
            # Strip the list of (for now) string values of their comma 
            dihedrals_list = list(map(lambda x: x.strip(","), dihedrals_list))
            if len(dihedrals_list) == 7:
                dihedrals_of_interest = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]
            else:
                print("Amount of dihedrals prompted is not aligned with the standard amount of dihedrals.")
                sys.exit(0)

            # Check if all values are float
            for i in dihedrals_list:
                if not isinstance(float(i), float):
                    print(f"One or more of the dihedral angles is not a floating point number : {i}. Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Check if size of the prompted dihedral values is the same as the amount of required dihedrals
            assert len(dihedrals_of_interest) == len(dihedrals_list),  "Check your input for missing dihedral values or missplaced commas.\nNote: the decimal values should be denoted by a point and not a comma."

            # Initialise dictionary
            set_of_dihedrals = {}
            # Append the values to their respective dihedrals 
            for dihr in range(len(dihedrals_of_interest)):
                set_of_dihedrals[dihedrals_of_interest[dihr]] = float(dihedrals_list[dihr])
            return set_of_dihedrals


        if moietyType == "linker":
            ## For now we hardcode this with the phospate linker, until we start broadening the linker space

            # Split the string into a list of strings
            dihedrals_list = list(map(lambda x : x.strip(","), dihedrals_list))

            # Check if all values are float
            for i in dihedrals_list:
                if not isinstance(float(i), float):
                    print(f"One or more of the dihedral angles is not a floating point number : {i}. Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Initialise dictionary
            set_of_dihedrals = {}
            set_of_dihedrals["dihedral_2"] = float(dihedrals_list[0])
            set_of_dihedrals["dihedral_1"] = float(dihedrals_list[1])
            return set_of_dihedrals


    def get_angles(self, identifier : str, moietyType : str, angles_list : list) -> dict:
        """ List the bond angles differently whether it belongs to a nucleoside or a linker moietyType """
        if moietyType == "nucleoside":
            angles_of_interest = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]

            # Strip the list of (for now) string values of their comma 
            angles_list = list(map(lambda x: x.strip(","), angles_list))

            # Check if all values are float
            for i in angles_list:
                if not isinstance(float(i), float):
                    print(f"One or more of the dihedral angles is not a floating point number : {i}. Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Check if size of the prompted dihedral values is the same as the amount of required dihedrals
            assert len(angles_of_interest) == len(angles_list),  "Check your input for missing dihedral values or missplaced commas.\nNote: the decimal values should be denoted by a point and not a comma."

            # Initialise dictionary
            set_of_angles = {}
            # Append the values to their respective bond angles
            for ang in range(len(angles_list)):
                set_of_angles[angles_of_interest[ang]] = float(angles_list[ang])
            return set_of_angles


        if moietyType == "linker":
            ## For now we hardcode this with the phospate linker, until we start broadening the linker space
            # Strip the list of (for now) string values of their comma 
            angles_list = list(map(lambda x: x.strip(","), angles_list))

            # Check if all values are float
            for i in angles_list:
                if not isinstance(float(i), float):
                    print(f"One or more of the dihedral angles is not a floating point number : {i}. Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Initialise dictionary
            set_of_angles = {}
            set_of_angles["OPO"] = float(angles_list[0])
            return set_of_angles


    def get_output_name(self, identifier : str, moietyType : str, conformation : Union[str, bool]) -> str:
        """ Create the name of the file based on the identifier of the nucleic acid chemistry and its corresponding base """
        if moietyType == "nucleoside":
            name_of_chemistry = identifier.lower()
            name_of_base = self.get_base().lower()

            if isinstance(conformation, bool):
                if not conformation:
                    return name_of_chemistry + "_" + name_of_base
            else :                                                      #is instance of string then
                conformation = conformation.lower()
                return name_of_chemistry + "_" + name_of_base + "_" + conformation

        elif moietyType == "linker":
            name_of_chemistry = identifier.lower()
            name_of_linker = TC.linker_dict[identifier].lower()
            return name_of_chemistry + "_" + name_of_linker

        else :
            sys.exit("The molecule is not annotated with either `nucleoside` or `linker`. Please revise the inputs")



class TransmuteToPdb:
    """ This class is used to convert the *.xyz files from ORCA to *.pdb files. Later on, these *.pdb files are prompted into Ducque to convert to *.json files.

        --xyz `*.xyz`
        --atomID `Residue Name`
        --atomname_list `see Ducque manual`               """


    def __init__(self, xyzfile):
        """ Initialise the object and create object attributes """
#        self.splitted = xyzfile.split('.')[0]
#        self.filename = self.splitted + '.xyz'
#        self.xyzname = self.basename + ".xyz"
#        self.array = np.array([])
#        self.pdb_dataframe = pd.DataFrame()
        self.pathName = xyzfile
        self.rootName = os.path.basename(xyzfile).split(".")[0]
        self.pdbName = self.rootName + ".pdb"
        self.x = list()
        self.y = list()
        self.z = list()
        self.elements = list()
        self.atomNameList = list()


    def parse_xyz_and_elementsymbol(self):
        """ Reads the xyz datafile and returns the coordinates and the element symbol"""
        # read only the lines with x-y-z coordinates
        with open(self.pathName, "r") as XYZ :
            _fileList = [line.strip() for line in XYZ.readlines()[2:]]
#            file = [line.strip() for line in file_list]
#        file_n = open(self.pathname, "r")
#        file_list = file_n.readlines()[2:]
#        file = [line.strip() for line in file_list]
#        file_n.close()

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
#        return x_coords, y_coords, z_coords, elements


    def return_processed_atomname_list(self, atomNameList : list):
        """ Returns a list of atom names in a well formatted list """
        atomNameList = list(map(lambda x : x.strip(","), atomNameList))

        for atom in atomNameList:
            if len(atom) > 4:
                print(f"The following atom has too many characters in the string {atom}. Maximum amount allowed is 4.\n")
                sys.exit(0)

        self.atomNameList = atomNameList


    def arraysize_vs_atomname_list_compatibility(self) -> bool:
        """ Checks to see if the size of the parsed cartesian coordinates (in length) matches the size of the inputted atomname_list. They should match! """

        if len(self.elements) == len(self.atomNameList):
            return True
        else :
            print(f"Atomname_list length : {len(self.atomNameList)}. Coordinate array size : ({len(self.elements)}, , 3 ).")
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
                            "Feel free to adjust this to your needs in the file : transmute_func.py ; elementsymbol_vs_atomname_list_compatibility() function.\n")
                    return False

            # First parse out only the first character of the atomNameList, since this the first letter denotes the element of that atom
            atom_type = self.atomNameList[i_atom][0]

            if atom_type != self.elements[i_atom]:
                print("Atoms that do not match : Atom prompted - {atom_type}. Element parsed : {elements[atom]}. Position {str(atom + 1)}\n")
                return False

        return True


    def write_to_pdb_format_file(self, atomID) :

        # Just create a random integer for the sequenceNumber, just needs to be filled in
        from random import randint
        randomised_integer_for_sequence_number = randint(1,100)


        atomNumbers = np.linspace(1, len(self.elements), len(self.elements), dtype=int)

        # Write out the *.pdb file
#        fileName = "testing_duplex.pdb"

        with open(self.pdbName, "w") as pdb:
            # Write out Leading Strand
            for idx in range(len(self.elements)):
                line = ["ATOM",
                        atomNumbers[idx],
                        self.atomNameList[idx],
                        atomID,
                        "A",
                        randomised_integer_for_sequence_number,
                        self.x[idx], self.y[idx],self.z[idx],
                        "1.00", "0.00",
                        self.elements[idx] ]

                pdb.write("%-4s  %5d %-4s %3s %s%4d    %8s%8s%8s%6s%6s          %2s\n" % tuple(line))



#    def fill_in_the_rest_of_the_pdb_dataframe_attribute(self, atomID, atomname_list, x_coords, y_coords, z_coords, elements):
#        """ Fill in the remaining blanks of the pdb dataframe"""
#        from random import randint
#        randomised_integer_for_sequence_number = randint(1,100)
#
#        AtomNum_range = np.linspace(1, len(elements), len(elements), dtype=int)
#
#        self.pdb_dataframe = pd.DataFrame(index=range(len(elements)))
#        self.pdb_dataframe['RecName'] = 'ATOM'
#        self.pdb_dataframe['AtomNum'] = AtomNum_range
#        self.pdb_dataframe['AtomName'] = atomname_list
#        self.pdb_dataframe['AltLoc'] = ' '
#        self.pdb_dataframe['ResName'] = atomID
#        self.pdb_dataframe['Chain'] = 'A'
#        self.pdb_dataframe['SeqNum'] = randomised_integer_for_sequence_number
#        self.pdb_dataframe['X_coord'] = x_coords
#        self.pdb_dataframe['Y_coord'] = y_coords
#        self.pdb_dataframe['Z_coord'] = z_coords
#        self.pdb_dataframe['Occupancy'] = '1.00'
#        self.pdb_dataframe['Temp'] = '0.00'
#        self.pdb_dataframe['SegmentID'] = '   '
#        self.pdb_dataframe['Element'] = elements
#
#
#    def write_to_pdb_formatted_file(self):
#        """ Write out the pdb file """
#
#        with open(self.splitted + ".pdb", "w") as pdb:
#        with open(self.pdbname, "w") as pdb:
#            for index, row in self.pdb_dataframe.iterrows():
#                split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
#                pdb.write("%-6s%5s%5s%s%3s%2s%5d  %8s%8s%9s%6s%7s%4s     %2s\n" % tuple(split_line))
#            pdb.write("END")
#            pdb.close()
