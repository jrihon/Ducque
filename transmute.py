import transmute_func as TF
import json


def Transmutation(pdb_file, nucleic_acid_chemistry : str, moiety : str, dihedral_list : list, angles_list : list):
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
    nucleic_acid = TF.TransmuteToJson(pdb_file)
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

    if moiety == "nucleoside":
        # full name
        fullname = nucleic_acid.get_full_name(nucleic_acid_chemistry, moiety)
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

    if moiety == "linker":
        # full name
        fullname = nucleic_acid.get_full_name(nucleic_acid_chemistry, moiety)
        identity.append(fullname)

        # abbreviated name
        abbr = nucleic_acid_chemistry
        identity.append(abbr)

    #------------------------------------- ANGLE ------------------------------------#
    # Initialise the dictionary for the dihedrals and the bond angles
    angle = {}

    # Get dihedrals
    dihedrals = nucleic_acid.get_dihedrals(nucleic_acid_chemistry, moiety, dihedral_list)

    # Get Bond angles
    angles = nucleic_acid.get_angles(nucleic_acid_chemistry, moiety, angles_list)

    angle["dihedrals"] = json.dumps(dihedrals)
    angle["bond_angles"] = json.dumps(angles)

    #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
    molecule["pdb_properties"] = pdb_properties
    molecule["identity"] = json.dumps(identity)
    molecule["angles"] = angle

    #----------------------------- WRITE OUT A JSON FILE ----------------------------#
    # The json dump() method always requires us to dump it to a file in the current directory
    fname = nucleic_acid.get_output_name(nucleic_acid_chemistry, moiety)

    with open("./json/" + fname + ".json", "w") as filejson:
        json.dump(molecule, filejson, indent=4)


