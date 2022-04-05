import sys
import json
from typing import Union

import sysDaedalus
import transmute_func as TF



def Transmutation(pdb_file, nucleicAcidChemistry : str, moiety : str, dihedralList : list, anglesList : list, conformation : Union[str, bool] = False ):
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
            The dihedrals of the backbone and the glycosidic dihedral
            The bond angles that correspond with the three last atoms of the respective dihedral

    """
    ## Read pdb and convert to dataframe
    nucleicAcid = TF.TransmuteToJson(pdb_file)
    nucleicAcid.pdb_for_attributes()


    ## Initialise the main dictionary
    molecule = {}

    #-------------------------------- PDB PROPERTIES --------------------------------#
    # Initialise pdb properties dictionary
    pdb_properties = {}

    # Coordinates and Shape
    pdb_properties["Coordinates"] = json.dumps(nucleicAcid.get_array())
    pdb_properties["Shape"] = json.dumps(nucleicAcid.get_shape_array())

    # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings, since there are apostrophes in the molecules
    pdb_properties["Atoms"] = json.dumps(nucleicAcid.get_atoms())

    # Get the element symbol
    pdb_properties["Symbol"] = json.dumps(nucleicAcid.get_element_symbol())

    #----------------------------------- IDENTITY -----------------------------------#
    # Initialise the identity list
    identity = []

    if moiety == "nucleoside":
        # full name
        _fullname = nucleicAcid.get_full_name(nucleicAcidChemistry, moiety)
        identity.append(_fullname)

        # abbreviated name
        _abbr = nucleicAcidChemistry
        identity.append(_abbr)

        # molecule chemistry
        _molecule_residuename = nucleicAcid.get_chemistry()
        identity.append(_molecule_residuename)

        # Base of the nucleic acid
        _base = nucleicAcid.get_base()
        identity.append(_base)

    if moiety == "linker":
        # full name
        _fullname = nucleicAcid.get_full_name(nucleicAcidChemistry, moiety)
        identity.append(_fullname)

        # abbreviated name
        _abbr = nucleicAcidChemistry
        identity.append(_abbr)

    #------------------------------ DIHEDRALS AND ANGLES -----------------------------#
    # Initialise the dictionary for the dihedrals and the bond angles
    angle = {}

    # Get dihedrals
    dihedrals = nucleicAcid.get_dihedrals(nucleicAcidChemistry, moiety, dihedralList)

    # Get Bond angles
    angles = nucleicAcid.get_angles(nucleicAcidChemistry, moiety, anglesList)

    angle["dihedrals"] = json.dumps(dihedrals)
    angle["bond_angles"] = json.dumps(angles)

    #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
    molecule["pdb_properties"] = pdb_properties
    molecule["identity"] = json.dumps(identity)
    molecule["angles"] = angle

    #----------------------------- WRITE OUT A JSON FILE ----------------------------#
    # The json dump() method always requires us to dump it to a file in the current directory
    _fname = nucleicAcid.get_output_name(nucleicAcidChemistry, moiety, conformation)

    # Get Daedalus home
    DAEDALUSHOME = sysDaedalus.return_DAEDALUS_home()

    # Write the inputfile to the Daedalus home directory
    with open(DAEDALUSHOME + "json/" + _fname + ".json", "w") as filejson:
        json.dump(molecule, filejson, indent=4)



def convert_XYZ_to_PDB(xyzFile, atomID : str, atomNameList : list):
    """ the main function that convert an xyz formatted file to the required pdb format """

    # Instantiate the object
    PdbToBe = TF.TransmuteToPdb(xyzFile)

    # Parse all the required data from the xyz file
    PdbToBe.parse_xyz_and_elementsymbol()
#    xCoords, yCoords, zCoords, elementSymbol = PDB_to_be.parse_xyz_and_elementsymbol()

    # Process the inputted atomname list from a string to a list
    PdbToBe.return_processed_atomname_list(atomNameList)
#    atomNameList = PDB_to_be.return_processed_atomname_list()

    # If the inputted atomname list and the array size do not match in size, exit the program
    if not PdbToBe.arraysize_vs_atomname_list_compatibility() :
        print("The size of the prompted '--atomname_list' is not equal to the array of the cartesian coordinates, pertaining to the atoms of the molecule.\n"
                "Please revise the prompted atomlist.")
        sys.exit(0)

    # If the inputted atomname list and the elements do not match in atoms, exit the program
    if not PdbToBe.elementsymbol_vs_atomname_list_compatibility() :
        print("The prompted '--atomname_list' do not match the order of the parsed ElementSymbol list, pertaining to the atoms of the molecule.\n"
                "Please revise the prompted atomlist.")
        sys.exit(0)

    # Write out the pdb file from all the gathered information
    PdbToBe.write_to_pdb_format_file(atomID)



#    PdbToBe.fill_in_the_rest_of_the_pdb_dataframe_attribute(atomID, atomNameList, xCoords, yCoords, zCoords, elementSymbol)
#
#    PdbToBe.write_to_pdb_formatted_file()
