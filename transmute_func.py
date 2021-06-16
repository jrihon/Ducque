import pandas as pd
import numpy as np
import json

import transmute_func_tools

class TransmuteToJson:

    def __init__(self, pdbfile):
        """ Initialise the object and create object properties"""
        self.splitted = pdbfile.split('.')[0]
        self.pdb_dataframe = pd.DataFrame()
        self.filename = self.splitted + '.pdb'
        self.array = np.array([])


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
        AtomName, ResName, Xcoord, Ycoord, Zcoord, ElementSym = ([] for i in range(6))

        # Read the file and fill out the dataframe
        with open(self.filename_readFile) as pdbfile:
            for line in pdbfile:
                if line[:4] == 'ATOM' or line[:6] == 'HETATM':

                    name = line[12:16]
                    AtomName.append(name)

                    resname = line[17:20]
                    ResName.append(resname)

                    _Xcoord = line[31:38]
                    Xcoord.append(_Xcoord)

                    _Ycoord = line[39:46]
                    Ycoord.append(_Ycoord)

                    _Zcoord = line[47:54]
                    Zcoord.append(_Zcoord)

                    ElemSym = line[76:78]
                    ElementSym.append(ElemSym)

            self.pdb_dataframe['AtomName'] = AtomName
            self.pdb_dataframe['ResName'] = ResName
            self.pdb_dataframe['X_Coord'] = Xcoord
            self.pdb_dataframe['Y_Coord'] = Ycoord
            self.pdb_dataframe['Z_Coord'] = Zcoord
            self.pdb_dataframe['ElementSymbol'] = ElementSym

            # Add the array as a property, to make it easier later on
            self.array = np.array([Xcoord, Ycoord, Zcoord], dtype=float).T


    def get_matrix(self) -> list:
        """ Create a list of the array of coordinates respectively."""
        return np.ndarray.tolist(self.array)


    def get_shape_array(self) -> list:
        """ Get the shape of the array, but in a list """
        return [self.array.shape[0], self.array.shape[1]]


    def get_atoms(self) -> list:
        """ Get atom names to a list, also strip any remaining whitespace in all the strings in the list."""
        atomlist = self.pdb_dataframe['AtomName'].tolist()
        return  list(map(lambda x: x.strip(), atomlist))


    def get_element_symbol(self):
        """ Get element symbol to a list, also strip any remaining whitespace in all the strings in the list."""
        elementlist = self.pdb_dataframe['ElementSymbol'].tolist()
        return  list(map(lambda x: x.strip(), elementlist))


    def get_ID(self) -> str:
        """ Get the name of the residue, also strip any remaining whitespace in all the strings in the list. """
        return self.pdb_dataframe['ResName'][0].strip()


    def get_full_name(identifier : str) -> str:
        """ Get the full name of the nucleic acid chemistry we want to convert to a json """
        return transmute_func_tools.identity_dict[identifier]


    def get_base(self) -> str:
        """ Get the base that corresponds with this nucleic acid. Take the last character of the string Residue Name and
            look for it in the dictionary """
        name_id = self.pdb_dataframe['ResName'][0]
        base = name_id[-1].upper()

        return transmute_func_tools.base_dict[base]


    def get_dihedrals(self, identifier : str) -> dict:
        # for alpha and zeta, we need to include the the atoms of the adjacent atoms
        # backbone_angles = ["P1", "O5'", "C5'", "C4'", "C3'", "O3'", "P2"]
        # we set the alpha and zeta angles fixed for now, since the atoms are \
                # not part of the nucleic acid residue
        dihedrals_of_interest = ["beta", "gamma", "delta", "epsilon", "chi"]

        # Initialise dictionary
        dihedral_dict = {}
        dihedral_dict["alpha"] = -39.202
        dihedral_dict["beta"] = pass
        dihedral_dict["gamma"] = pass
        dihedral_dict["delta"] = pass
        dihedral_dict["epsilon"] = pass
        dihedral_dict["zeta"] = -98.887

        dihedral_dict["chi"] = pass

        return dihedral_dict

    def get_angles(self, identifier : str) -> dict:

        # Initialise dictionary
        angle_dict = {}


        return angle_dict


    def get_file_name(self, identifier : str) -> str:
        """ Create the name of the file based on the identifier of the nucleic acid chemistry and its corresponding base """
        name_of_base = self.get_base()

##---------------------------- FUNCTIONS THAT ARE NOT IN USE ANYMORE ----------------------------##
#class TransmuteToPDB:
#
#    def __init__(self, jsonfile):
#
#        self.splitted = jsonfile.split('.')[0]
#        self.pdb_dataframe = pd.DataFrame()
#        self.filename = self.splitted + '.json'
#
#        with open(self.filename) as jason:
#            self.jason = json.load(jason)
#
#
#    def get_recordName(self):
#
#        length_array = json.loads(self.jason['pdb_properties']['Shape'])[0]
#        return ['ATOM' for x in range(length_array)]
#
#
#    def get_coord_array(self):
#
#        # Since the dimensions are preserved, we don't need to worry about the shape of the array
#        return np.asarray(json.loads(self.jason['pdb_properties']['Coordinates']), dtype=float)
#
#
#    def get_ID(self):
#
#        return json.loads(self.jason['Identity'])[2]
#
#
#    def get_atoms(self):
#
#        return json.loads(self.jason['pdb_properties']['Atoms'])
#
#
#    def get_sequence(self):
#
#        length_sequence = json.loads(self.jason['pdb_properties']['Shape'])[0]
#        return [x for x in range(1, length_sequence + 1)]
#
#
#    def get_symbol(self):
#
#        return json.loads(self.jason['pdb_properties']['Symbol'])
#
#
#    def write_pdb(self):
#        with open(self.splitted + '.pdb' ,'w') as pdb:
#            for index, row in self.pdb_dataframe.iterrows():
#                split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
#                pdb.write('%-6s%5s%5s%s%s%2s%5s   %8s%8s%8s%6s%6s%4s      %2s\n' % tuple(split_line))
#
#            pdb.write('END')
#
### just \n ##
