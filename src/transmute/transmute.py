import sys
import json
from typing import Union

import systemsDucque
import transmute.utils_transmute as UT



def Transmutation(pdb_fname, nucleicAcidChemistry : str, moietyType : str, dihedralList : list, anglesList : list, conformation : Union[str, bool] = False ):
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
    nucleicAcid = UT.TransmuteToJson(pdb_fname)
    nucleicAcid.pdb_for_attributes()


    ## Initialise the main dictionary
    molecule = {}

    #-------------------------------- PDB PROPERTIES --------------------------------#
    # Initialise pdb properties dictionary
    pdb_properties = {}

    # Coordinates and Shape
    pdb_properties["Coordinates"] = json.dumps(nucleicAcid.get_array())
    pdb_properties["Shape"] = json.dumps(nucleicAcid.get_shape_array())

    # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings (escape characters), since there are apostrophes in the molecules
    pdb_properties["Atoms"] = json.dumps(nucleicAcid.get_atoms())

    # Get the element symbol
    pdb_properties["Symbol"] = json.dumps(nucleicAcid.get_element_symbol())

    #----------------------------------- IDENTITY -----------------------------------#
    # Initialise the identity list
    identity = []

    if moietyType == "nucleoside":
        # full name
        _fullname = nucleicAcid.get_full_name(nucleicAcidChemistry, moietyType)
        identity.append(_fullname)

        # abbreviated name
        _abbr = nucleicAcidChemistry
        identity.append(_abbr)

        # molecule chemistry, which is often the same as the residue name in the pdb
        _molecule_residuename = nucleicAcid.get_chemistry()
        identity.append(_molecule_residuename)

        # Nucleobase of the nucleic acid
        _nucleobase = nucleicAcid.get_base()
        identity.append(_nucleobase)

    if moietyType == "linker":
        # full name
        _fullname = nucleicAcid.get_full_name(nucleicAcidChemistry, moietyType)
        identity.append(_fullname)

        # abbreviated name
        _abbr = nucleicAcidChemistry
        identity.append(_abbr)

    #------------------------------- TORSIONS AND ANGLES -----------------------------#
    # Initialise the dictionary for the dihedrals and the bond angles
    angles = {}

    # Get dihedrals
    _torsions = nucleicAcid.get_dihedrals(nucleicAcidChemistry, moietyType, dihedralList)

    # Get Bond angles
    _angles = nucleicAcid.get_angles(nucleicAcidChemistry, moietyType, anglesList)

    angles["dihedrals"] = json.dumps(_torsions)
    angles["bond_angles"] = json.dumps(_angles)

    #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
    molecule["pdb_properties"] = pdb_properties
    molecule["identity"] = json.dumps(identity)
    molecule["angles"] = angles

    #----------------------------- WRITE OUT A JSON FILE ----------------------------#
    # The json dump() method always requires us to dump it to a file in the current directory
    _fname = nucleicAcid.get_output_name(nucleicAcidChemistry, moietyType, conformation)

    # Get Ducque home
    DUCQUEHOME = systemsDucque.return_DUCQUEHOME()

    # Write the inputfile to the Ducque home directory
    with open(DUCQUEHOME + "json/" + _fname + ".json", "w") as filejson:
        json.dump(molecule, filejson, indent=4)



def convert_XYZ_to_PDB(xyzFname : str, atomID : str, atomNameList : list):
    """ the main function that convert an xyz formatted file to the required pdb format """

    # Instantiate the object
    PdbToBe = UT.TransmuteToPdb(xyzFname)

    # Parse all the required data from the xyz file
    PdbToBe.parse_xyz_and_elementsymbol()
#    xCoords, yCoords, zCoords, elementSymbol = PDB_to_be.parse_xyz_and_elementsymbol()

    # Process the inputted atomname list from a string to a list
    PdbToBe.return_processed_atomname_list(atomNameList)

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
