import transmute_func
import json


def NA_to_json(pdb_file):
    ## Initialise the dictionary
    molecule = {}
    pdb_properties = {}
    # backbone is going to be a part of molecule; nested
    backbone = {}
    Identity = ['Deoxyribonucleic acid', 'DNA']

    ## Read pdb and convert to dataframe
    nucleic_acid = transmute_func.TransmuteToJson(pdb_file)
    nucleic_acid.pdb_to_dataframe()

    ## Retrieve matrix about the coordinates
    xyz_array, shape = nucleic_acid.write_matrix()
    pdb_properties['Coordinates'], pdb_properties['Shape'] = json.dumps(xyz_array), json.dumps(shape)

    ## Get Atom namelist
    # NB: json outputs double quotations as \" XX \" for strings, \
            # since there are apostrophes in the molecules
    pdb_properties['Atoms'] = json.dumps(nucleic_acid.write_atoms())

    ## Get the dihedral angles
    # Here we'll ask a list to be prompted, that contain the backbone atoms in the correct order.
    #backbone_atoms = ["P1", "O5'", "C5'", "C4'", "C3'", "O3'", "P2"]
    backbone_dihr = nucleic_acid.write_backbone_dihedrals()

    backbone['Backbone'] = json.dumps(backbone_dihr)

    # Element symbol
    pdb_properties['Symbol'] = json.dumps(nucleic_acid.write_element_symbol())

    # molecule ID
    Identity.append(nucleic_acid.write_ID())

    ## Put everything into the pdb
    molecule['pdb_properties'] = pdb_properties
    molecule['Identity'] = json.dumps(Identity)
    molecule['Dihedrals'] = backbone

    # The json dump() method always requires us to dump it to a file in the current directory
    fname = nucleic_acid.filename.split('/')[-1].split('.')[0]
    filejson = open(fname + '.json', 'w')
    dump_json = json.dump(molecule, filejson, indent=4)
    filejson.close()


def NA_to_pdb(json_file):

    # Create object
    nucleic_acid = transmute_func.TransmuteToPDB(json_file)

    # See if it is possible to extracts the array correctly
    coord_array = nucleic_acid.get_coord_array()

    # Create dataframe
    nucleic_acid.pdb_dataframe['RecName'] = nucleic_acid.get_recordName()
    nucleic_acid.pdb_dataframe['AtomNum'] = nucleic_acid.get_sequence()
    nucleic_acid.pdb_dataframe['AtomName'] = nucleic_acid.get_atoms()
    nucleic_acid.pdb_dataframe['AltLoc'] = ' '
    nucleic_acid.pdb_dataframe['ResName'] = nucleic_acid.get_ID()
    nucleic_acid.pdb_dataframe['Chain'] = 'A'
    nucleic_acid.pdb_dataframe['Sequence'] = str(1)
    nucleic_acid.pdb_dataframe['X_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,0]))
    nucleic_acid.pdb_dataframe['Y_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,1]))
    nucleic_acid.pdb_dataframe['Z_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,2]))
    nucleic_acid.pdb_dataframe['Occupancy'] = '1.00'
    nucleic_acid.pdb_dataframe['Temp'] = '0.00'
    nucleic_acid.pdb_dataframe['SegmentID'] = str('   ')
    nucleic_acid.pdb_dataframe['ElementSym'] = nucleic_acid.get_symbol()

    nucleic_acid.write_pdb()


def linker_to_json(pdb_file):

    ## Initialise the dictionary
    molecule = {}
    pdb_properties = {}
    # dihedrals of the phosphate ( C5', O5' , P , OP1)
    dihedrals = {'dihedral_oxygen_OP1': 76.676, 'dihedral_oxygen_OP2': -155.039}
    angles = {'C5_O5_P' : 118.958, 'O5_P_OP' : 109.766, 'O5_P_O3' : 101.415}
    Identity = ['Phosphate','Linker']

    # Create instance
    linker = transmute_func.TransmuteToJson(pdb_file)

    linker.pdb_to_dataframe()

    # Retrieve matrix about coordinates
    xyz_array, shape = linker.write_matrix()
    pdb_properties['Coordinates'], pdb_properties['Shape'] = json.dumps(xyz_array), json.dumps(shape)

    # Atom namelist
    pdb_properties['Atoms'] = json.dumps(linker.write_atoms())

    # Element Symbol
    pdb_properties['Symbol'] = json.dumps(linker.write_element_symbol())

    # Put everything into the pdb
    molecule['pdb_properties'] = pdb_properties
    molecule['Identity'] = Identity
    molecule['Dihedrals'] = dihedrals
    molecule['Angles'] = angles


    # The json dump() method always requires us to dump it to a file in the current directory
    fname = linker.filename.split('/')[-1].split('.')[0]
    filejson = open(fname + '.json', 'w')
    dump_json = json.dump(molecule, filejson, indent=4)
    filejson.close()


def linker_to_pdb(json_file):

    # Create object
    linker = transmute_func.TransmuteToPDB(json_file)

    # See if it is possible to extracts the array correctly
    coord_array = linker.get_coord_array()

    # Create dataframe
    linker.pdb_dataframe['RecName'] = linker.get_recordName()
    linker.pdb_dataframe['AtomNum'] = linker.get_sequence()
    linker.pdb_dataframe['AtomName'] = linker.get_atoms()
    linker.pdb_dataframe['AltLoc'] = ' '
    linker.pdb_dataframe['ResName'] = 'POS'     # This is temporary, just to see what the results are
    linker.pdb_dataframe['Chain'] = 'A'
    linker.pdb_dataframe['Sequence'] = str(1)
    linker.pdb_dataframe['X_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,0]))
    linker.pdb_dataframe['Y_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,1]))
    linker.pdb_dataframe['Z_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,2]))
    linker.pdb_dataframe['Occupancy'] = '1.00'
    linker.pdb_dataframe['Temp'] = '0.00'
    linker.pdb_dataframe['SegmentID'] = str('   ')
    linker.pdb_dataframe['ElementSym'] = linker.get_symbol()

    linker.write_pdb()


def pdb_is_linker(pdb_file):
    """ Check to see if the pdb file that is loaded in is a linker molecule"""

    with open(pdb_file) as pdb:
        lines = len(pdb.readlines())

    # if the amount of lines is presumably larger than the amount of atoms in a linker segment
    if lines > 10:
        return False
    else:
        return True


def json_is_linker(json_file):
    """ Check to see if the json file that is loaded in is a linker molecule """

    # Load json data into the instance 'check'
    check = transmute_func.TransmuteToPDB(json_file)

    if check.jason['Identity'][1] == 'Linker':
        return True
    else:
        return False






















## just \n ##
