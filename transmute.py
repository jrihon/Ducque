import transmute_func as TF
import json


def Transmutation(pdb_file, nucleic_acid_chemistry : str):
    """This function converts a pdb formatted file into a json file.
    Json files make for a much easier data parsing format, are computationally much more efficient and require less memory to be held.

    Json Format:
        molecule:
            The properties parsed from the pdb_file
            The shape of the coordinate array
            The name of the respective atoms of the molecule
            The element symbol of the respective atoms of the molecule
        identity:
            [the name of the chemistry, abbreviated name of the chemistry, residue name, type of base]
        Angle:
            The dihedrals of the backbone and the anomeric carbon-base
            The bond angles that correspond with the three last atoms of the respective dihedral

    """
    ## Read pdb and convert to dataframe
    nucleic_acid = transmute_func.TransmuteToJson(pdb_file)
    nucleic_acid.pdb_to_dataframe()


    ## Initialise the main dictionary
    molecule = {}

    #-------------------------------- PDB PROPERTIES --------------------------------#
    # Initialise pdb properties dictionary
    pdb_properties = {}

    # Coordinates and Shape
    pdb_properties["Coordinates"] = json.dumps(nucleic_acid.get_array())
    pdb_properties["Shape"] = json.dumps(nucleic_acid.get_shape_array())

    # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings, since there are apostrophes in the molecules
    pdb_properties["Atoms"] = json.dumps(nucleic_acid.get_atoms())

    # Get the element symbol
    pdb_properties["Symbol"] = json.dumps(nucleic_acid.get_element_symbol())

    #----------------------------------- IDENTITY -----------------------------------#
    # Initialise the identity list
    identity = []

    # full name
    fullname = TF.get_full_name(nucleic_acid_chemistry)
    identity.append(fullname)

    # abbreviated name
    abbr = nucleic_acid_chemistry
    identity.append(abbr)

    # molecule ID
    molecule_residuename = nucleic_acid.get_ID()
    identity.append(molecule_residuename)

    # Base of the nucleic acid
    base = nucleic_acid.get_base()
    identity.append(base)

    #------------------------------------- ANGLE ------------------------------------#
    # Initialise the dictionary for the dihedrals and the bond angles
    angle = {}

    # Get dihedrals
    dihedrals = nucleic_acid.get_dihedrals(nucleic_acid_chemistry, base)

    # Get Bond angles
    angles = nucleic_acid.get_angles(nucleic_acid_chemistry, base)

    angle["dihedrals"] = json.dumps(dihedrals)
    angle["bond_angles"] = json.dumps(angles)

    #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
    molecule["pdb_properties"] = pdb_properties
    molecule["identity"] = json.dumps(identity)
    molecule["angles"] = angle

    #----------------------------- WRITE OUT A JSON FILE ----------------------------#
#    # The json dump() method always requires us to dump it to a file in the current directory
#    fname = nucleic_acid.get_file_name(nucleic_acid_chemistry)
#
#    with open("./json/" + fname + ".json", "w"):
#        dump_json = json.dump(molecule, filejson, indent=4)


def linker_to_json(pdb_file):

    ## Initialise the dictionary
    molecule = {}
    pdb_properties = {}
    # dihedrals of the phosphate ( C5', O5' , P , OP1)
    dihedrals = {'dihedral_oxygen_OP1': 76.676, 'dihedral_oxygen_OP2': -155.039}
    angles = {'C5_O5_P' : 118.958, 'O5_P_OP2' : 109.766, 'O5_P_OP1' : 109.746, 'O5_P_O3' : 101.415}
    identity = ['Phosphate','Linker']

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
    molecule['identity'] = identity
    molecule['Dihedrals'] = dihedrals
    molecule['Angles'] = angles


    # The json dump() method always requires us to dump it to a file in the current directory
    fname = linker.filename.split('/')[-1].split('.')[0]
    filejson = open(fname + '.json', 'w')
    dump_json = json.dump(molecule, filejson, indent=4)
    filejson.close()


##---------------------------- FUNCTIONS THAT ARE NOT IN USE ANYMORE ----------------------------###
#def NA_to_pdb(json_file):
#    # Create object
#    linker = transmute_func.TransmuteToPDB(json_file)
#
#    # See if it is possible to extracts the array correctly
#    coord_array = linker.get_coord_array()
#
#    # Create dataframe
#    linker.pdb_dataframe['RecName'] = linker.get_recordName()
#    linker.pdb_dataframe['AtomNum'] = linker.get_sequence()
#    linker.pdb_dataframe['AtomName'] = linker.get_atoms()
#    linker.pdb_dataframe['AltLoc'] = ' '
#    linker.pdb_dataframe['ResName'] = 'POS'     # This is temporary, just to see what the results are
#    linker.pdb_dataframe['Chain'] = 'A'
#    linker.pdb_dataframe['Sequence'] = str(1)
#    linker.pdb_dataframe['X_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,0]))
#    linker.pdb_dataframe['Y_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,1]))
#    linker.pdb_dataframe['Z_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,2]))
#    linker.pdb_dataframe['Occupancy'] = '1.00'
#    linker.pdb_dataframe['Temp'] = '0.00'
#    linker.pdb_dataframe['SegmentID'] = str('   ')
#    linker.pdb_dataframe['ElementSym'] = linker.get_symbol()
#
#    linker.write_pdb()


#def linker_to_pdb(json_file):
#
#    # Create object
#    linker = transmute_func.TransmuteToPDB(json_file)
#
#    # See if it is possible to extracts the array correctly
#    coord_array = linker.get_coord_array()
#
#    # Create dataframe
#    linker.pdb_dataframe['RecName'] = linker.get_recordName()
#    linker.pdb_dataframe['AtomNum'] = linker.get_sequence()
#    linker.pdb_dataframe['AtomName'] = linker.get_atoms()
#    linker.pdb_dataframe['AltLoc'] = ' '
#    linker.pdb_dataframe['ResName'] = 'POS'     # This is temporary, just to see what the results are
#    linker.pdb_dataframe['Chain'] = 'A'
#    linker.pdb_dataframe['Sequence'] = str(1)
#    linker.pdb_dataframe['X_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,0]))
#    linker.pdb_dataframe['Y_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,1]))
#    linker.pdb_dataframe['Z_coord'] = list(map(lambda x: '{:.3f}'.format(x), coord_array[:,2]))
#    linker.pdb_dataframe['Occupancy'] = '1.00'
#    linker.pdb_dataframe['Temp'] = '0.00'
#    linker.pdb_dataframe['SegmentID'] = str('   ')
#    linker.pdb_dataframe['ElementSym'] = linker.get_symbol()
#
#    linker.write_pdb()
#
#def pdb_is_linker(pdb_file):
#    """ Check to see if the pdb file that is loaded in is a linker molecule"""
#
#    with open(pdb_file) as pdb:
#        lines = len(pdb.readlines())
#
#    # if the amount of lines is presumably larger than the amount of atoms in a linker segment
#    if lines > 10:
#        return False
#    else:
#        return True
#
#
#def json_is_linker(json_file):
#    """ Check to see if the json file that is loaded in is a linker molecule """
#
#    # Load json data into the instance 'check'
#    check = transmute_func.TransmuteToPDB(json_file)
#
#    if check.jason['identity'][1] == 'Linker':
#        return True
#    else:
#        return False

## just \n ##
