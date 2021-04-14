import pandas as pd
import numpy as np
#import pytraj as ptj
import json


class TransmuteToJson:

    def __init__(self, pdbfile):
        """
        Initialise a dataframe 
        """
        self.splitted = pdbfile.split('.')[0]
        self.pdb_dataframe = pd.DataFrame()
        self.filename = self.splitted + '.pdb'

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
        # Initialize lists as a 
        AtomName, ResName, Xcoord, Ycoord, Zcoord, ElementSym = ([] for i in range(6))

        _readFile = self.splitted + '.pdb'
        with open(_readFile) as pdbfile:
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


    def write_matrix(self):

        # Get array
        Xarray = self.pdb_dataframe['X_Coord']
        Yarray = self.pdb_dataframe['Y_Coord']
        Zarray = self.pdb_dataframe['Z_Coord']

        # Transpose the arraym so it can get properly shaped
        coord_array = np.array([Xarray, Yarray, Zarray], dtype=float).T
        coord_list = np.ndarray.tolist(coord_array)

        # Get array shape
        array_shape = coord_array.shape

        # Return the values 
        return coord_list, array_shape


    def write_ID(self):
        # Get the name of the residue
        return self.pdb_dataframe['ResName'][0].strip()


    def write_atoms(self):

        # Get it to a list
        atomlist = self.pdb_dataframe['AtomName'].tolist()

        return  list(map(lambda x: x.strip(), atomlist))


    def write_backbone_dihedrals(self):
        # for alpha and zeta, we need to include the the atoms of the adjacent atoms
        # backbone_angles = ["P1", "O5'", "C5'", "C4'", "C3'", "O3'", "P2"]
        # we set the alpha and zeta angles fixed for now, since the atoms are \
                # not part of the nucleic acid residue
        backbone_dihedrals = ['beta', 'gamma', 'delta', 'epsilon']

        # Load the molecule
#        mol_nucleicacid = ptj.load(self.filename)
        # Initialise dictionary
        backbone_dict = {}
        backbone_dict['alpha'] = -39.202
#        for i in range(len(backbone_dihedrals)):
#            # Create atom mask suitable for pytraj to read it in
#            atom_mask = list(map(lambda x: '@' + x, [ BA[i], BA[i + 1], BA[i + 2], BA[i + 3] ]))
#            atom_mask = ', '.join(atom_mask)
#
#            # Calculate dihedral angle and extract and append immediately to the dictionary
#            _key = backbone_dihedrals[i]
#            _dihedr = ptj.dihedral(traj=mol_nucleicacid, mask=atom_mask)
#            
#            backbone_dict[_key] = round(_dihedr[0], 3)
#
#
        backbone_dict['beta'] = -151.431
        backbone_dict['gamma'] = 30.929
        backbone_dict['delta'] = 156.517
        backbone_dict['epsilon'] = 159.171

        backbone_dict['zeta'] = -98.887
        return backbone_dict


    def write_element_symbol(self):
        
        # Get it to a list 
        atomlist = self.pdb_dataframe['ElementSymbol'].tolist()
        return  list(map(lambda x: x.strip(), atomlist))



class TransmuteToPDB:

    def __init__(self, jsonfile):

        self.splitted = jsonfile.split('.')[0]
        self.pdb_dataframe = pd.DataFrame()
        self.filename = self.splitted + '.json'

        with open(self.filename) as jason:
            self.jason = json.load(jason)


    def get_recordName(self):

        length_array = json.loads(self.jason['pdb_properties']['Shape'])[0]
        return ['ATOM' for x in range(length_array)]


    def get_coord_array(self):

        # Since the dimensions are preserved, we don't need to worry about the shape of the array
        return np.asarray(json.loads(self.jason['pdb_properties']['Coordinates']), dtype=float)


    def get_ID(self):

        return json.loads(self.jason['Identity'])[2]


    def get_atoms(self):

        return json.loads(self.jason['pdb_properties']['Atoms'])


    def get_sequence(self):

        length_sequence = json.loads(self.jason['pdb_properties']['Shape'])[0]
        return [x for x in range(1, length_sequence + 1)]


    def get_symbol(self):

        return json.loads(self.jason['pdb_properties']['Symbol'])


    def write_pdb(self):
        with open(self.splitted + '.pdb' ,'w') as pdb:
            for index, row in self.pdb_dataframe.iterrows():
                split_line = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13] ]
                pdb.write('%-6s%5s%5s%s%s%2s%5s   %8s%8s%8s%6s%6s%4s      %2s\n' % tuple(split_line))

            pdb.write('END')

## just \n ##
