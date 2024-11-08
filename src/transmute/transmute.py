import json
from os import getcwd, listdir
from os.path import isfile, isdir

import systemsDucque as SD
import transmute.utils_transmute as UT



def Transmutation(pdb_fname, nucleicAcidChemistry : str, moietyType : str, dihedralList : list, anglesList : list, conformation : str , nucleobase : str):
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

    # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings (escape characters)
    pdb_properties["Atoms"] = json.dumps(nucleicAcid.get_atoms())

    # Get the element symbol
    pdb_properties["Symbol"] = json.dumps(nucleicAcid.get_element_symbol())

    #----------------------------------- IDENTITY -----------------------------------#
    # Initialise the identity list
    identity = []

    if moietyType.upper() == "NUCLEOSIDE":
        # full name
        fullname = nucleicAcid.get_full_name(nucleicAcidChemistry, moietyType)
        identity.append(fullname)

        # abbreviated name
#        abbr = nucleicAcidChemistry
        identity.append(nucleicAcidChemistry)

        # molecule chemistry, which is often the same as the residue name in the pdb
        molecule_residuename = nucleicAcid.get_resname()
        identity.append(molecule_residuename)

        # Nucleobase of the nucleic acid
        nucleobaseName = nucleicAcid.get_nucleobase(nucleobase)
        identity.append(nucleobaseName)

        # Check if atoms to build by, stated in the TABLES, are present in the pdb
        nucleicAcid.validate_atomnames_for_building(moietyType=moietyType, chemistry=nucleicAcidChemistry, nucleobase=nucleobaseName)

    if moietyType.upper() == "LINKER":
        # full name
        fullname = nucleicAcid.get_full_name(nucleicAcidChemistry, moietyType)
        identity.append(fullname)

        # abbreviated name
#        abbr = nucleicAcidChemistry
        identity.append(nucleicAcidChemistry)

        # Check if atoms to build by, stated in the TABLES, are present in the pdb
        print(nucleicAcidChemistry)
        nucleicAcid.validate_atomnames_for_building(moietyType=moietyType, chemistry=nucleicAcidChemistry)

    #------------------------------- TORSIONS AND ANGLES -----------------------------#
    # Initialise the dictionary for the dihedrals and the bond angles
    angles = {}

    # Get Bond angles
    bondangles = nucleicAcid.get_angles(moietyType, anglesList)

    # Get dihedrals
    torsions = nucleicAcid.get_angles( moietyType, dihedralList)

    angles["bond_angles"] = json.dumps(bondangles)
    angles["dihedrals"] = json.dumps(torsions)

    #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
    molecule["pdb_properties"] = pdb_properties
    molecule["identity"] = json.dumps(identity)
    molecule["angles"] = angles

    #----------------------------- WRITE OUT A JSON FILE ----------------------------#
    # The json dump() method always requires us to dump it to a file in the current directory
    fname = nucleicAcid.get_output_name(nucleicAcidChemistry, moietyType, conformation, nucleobase)

    # Get Ducque home
    DUCQUEHOME = SD.return_DUCQUEHOME()

    # Write the inputfile to the Ducque home directory
    with open(DUCQUEHOME + "json/" + fname + ".json", "w") as filejson:
        json.dump(molecule, filejson, indent=4)

    SD.print_writing(f"{DUCQUEHOME}json/{fname}.json")



def convert_XYZ_to_PDB(xyzFname : str, residue : str, atomNameList : list):
    """ the main function that convert an xyz formatted file to the required pdb format """

    pathToFile = getcwd() + "/" + xyzFname

    # If is file
    if isfile(pathToFile) :
        # Instantiate the object
        PdbToBe = UT.TransmuteToPdb(xyzFname)

        # Parse all the required data from the xyz file
        PdbToBe.parse_xyz_and_elementsymbol()

        # Process the inputted atomname list from a string to a list
        PdbToBe.return_processed_atomname_list(atomNameList)

        # If the inputted atomname list and the array size do not match in size, exit the program
        if not PdbToBe.arraysize_vs_atomname_list_compatibility() :
            print("The size of the prompted '--atomname_list' is not equal to the array of the cartesian coordinates, pertaining to the atoms of the molecule.\n"
                    "Please revise the prompted atomlist.")
            SD.exit_Ducque()

        # If the inputted atomname list and the elements do not match in atoms, exit the program
        if not PdbToBe.elementsymbol_vs_atomname_list_compatibility() :
            print("The prompted '--atomname_list' do not match the order of the parsed ElementSymbol list, pertaining to the atoms of the molecule.\n"
                    "Please revise the prompted atomlist.")
            SD.exit_Ducque()

        # Write out the pdb file from all the gathered information
        PdbToBe.write_to_pdb_format_file(residue)

    # If is directory
    elif isdir(pathToFile):
        # find all xyz files in a directory
        lsdir = [x for x in listdir(pathToFile) if x.endswith(".xyz") and not x.endswith("_trj.xyz")]  # xyz coordinate files, not FNAME_trj.xyz

        # Do a first check
        PdbToBe = UT.TransmuteToPdb(xyzFname + "/" + lsdir[0])
        PdbToBe.parse_xyz_and_elementsymbol()
        PdbToBe.return_processed_atomname_list(atomNameList)

        if not PdbToBe.arraysize_vs_atomname_list_compatibility() :
            print("The size of the prompted '--atomname_list' is not equal to the array of the cartesian coordinates, pertaining to the atoms of the molecule.\n"
                    "Please revise the prompted atomlist.")
            SD.exit_Ducque()

        if not PdbToBe.elementsymbol_vs_atomname_list_compatibility() :
            print("The prompted '--atomname_list' do not match the order of the parsed ElementSymbol list, pertaining to the atoms of the molecule.\n"
                    "Please revise the prompted atomlist.")
            SD.exit_Ducque()

        counter = 0
        # If this all works, don't do the compatibility checks anymore
        for fname in lsdir: 
            PdbToBe = UT.TransmuteToPdb(xyzFname + "/" + fname)
            PdbToBe.parse_xyz_and_elementsymbol()
            PdbToBe.return_processed_atomname_list(atomNameList)
            PdbToBe.write_to_pdb_format_file(residue)
            counter += 1

        SD.print_writing(f"Written to {counter} pdb formatted files")

    else : 
        SD.print_filenotfound(pathToFile)
        
