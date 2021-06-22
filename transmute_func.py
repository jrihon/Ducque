import pandas as pd
import numpy as np
import json, sys, os

import transmute_func_tools as TFT

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

        # Check if file is in cwd or in the pdb directory
        pdbfname = self.filename
        if not os.path.isfile(pdbfname):
            pdbfname = "./pdb/" + self.filename

        # Read the file and fill out the dataframe
        with open(pdbfname) as pdbfile:
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

            # Add the array as an attribute
            self.array = np.array([Xcoord, Ycoord, Zcoord], dtype=float).T

            # Add the atom name list as an attribute
            self.atom_list = list(map(lambda x : x.strip(), AtomName))

    def get_array(self) -> list:
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


    def get_full_name(self, identifier : str, moiety) -> str:
        """ Get the full name of the nucleic acid chemistry we want to convert to a json """
        if moiety == "nucleoside":
            return TFT.nucleoside_dict[identifier.upper()]
        if moiety == "linker":
            return TFT.linker_dict[identifier.upper()]


    def get_base(self) -> str:
        """ Get the base that corresponds with this nucleic acid. Take the last character of the string Residue Name and
            look for it in the dictionary """
        name_id = self.pdb_dataframe['ResName'][0].strip()
        base = str(name_id)[-1]

        return TFT.base_dict[base]


    def get_dihedrals(self, identifier : str, moiety : str,  dihedrals_list : list) -> dict:
        """ List the dihedrals differently whether it belongs to a nucleoside or a linker moiety """

        if moiety == "nucleoside":
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
                    print("One or more of the dihedral angles is not a floating point number : " + i + ". Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Initialise dictionary
            set_of_dihedrals = {}
            # Append the values to their respective dihedrals 
            for dihr in range(len(dihedrals_of_interest)):
                set_of_dihedrals[dihedrals_of_interest[dihr]] = float(dihedrals_list[dihr])
            return set_of_dihedrals


        if moiety == "linker":
            ## For now we hardcode this with the phospate linker, until we start broadening the linker space

            # Split the string into a list of strings
            dihedrals_list = list(map(lambda x : x.strip(","), dihedrals_list))

            # Check if all values are float
            for i in dihedrals_list:
                if not isinstance(float(i), float):
                    print("One or more of the dihedral angles is not a floating point number : " + i + ". Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Initialise dictionary
            set_of_dihedrals = {}
            set_of_dihedrals["OP2_dihedral"] = float(dihedrals_list[0])
            set_of_dihedrals["OP1_dihedral"] = float(dihedrals_list[1])
            return set_of_dihedrals


    def get_angles(self, identifier : str, moiety : str, angles_list : list) -> dict:
        """ List the bond angles differently whether it belongs to a nucleoside or a linker moiety """
        if moiety == "nucleoside":
            angles_of_interest = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]

            # Strip the list of (for now) string values of their comma 
            angles_list = list(map(lambda x: x.strip(","), angles_list))

            # Check if all values are float
            for i in angles_list:
                if not isinstance(float(i), float):
                    print("One or more of the dihedral angles is not a floating point number : " + i + ". Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Initialise dictionary
            set_of_angles = {}
            # Append the values to their respective bond angles
            for ang in range(len(angles_list)):
                set_of_angles[angles_of_interest[ang]] = float(angles_list[ang])
            return set_of_angles


        if moiety == "linker":
            ## For now we hardcode this with the phospate linker, until we start broadening the linker space
            # Strip the list of (for now) string values of their comma 
            angles_list = list(map(lambda x: x.strip(","), angles_list))

            # Check if all values are float
            for i in angles_list:
                if not isinstance(float(i), float):
                    print("One or more of the dihedral angles is not a floating point number : " + i + ". Please reconsider the entries for the dihedrals.\n")
                    sys.exit(0)

            # Initialise dictionary
            set_of_angles = {}
            set_of_angles["OPO"] = float(angles_list[0])
            return set_of_angles


    def get_output_name(self, identifier : str, moiety : str) -> str:
        """ Create the name of the file based on the identifier of the nucleic acid chemistry and its corresponding base """
        if moiety == "nucleoside":
            name_of_chemistry = identifier.lower()
            name_of_base = self.get_base().lower()
            return name_of_chemistry + "_" + name_of_base

        if moiety == "linker":
            name_of_chemistry = identifier.lower()
            name_of_linker = TFT.linker_dict[identifier].lower()
            return name_of_chemistry + "_" + name_of_linker

