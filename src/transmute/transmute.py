import json
from os import getcwd, listdir
from os.path import isfile, isdir

import systemsDucque as SD
import transmute.utils_transmute as UT



def Transmutation(transmuteObject : UT.NucleosideTransmute | UT.LinkerTransmute | UT.NucleobaseTransmute) -> None :
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

        Admittedly quite a large function, but I'll take readability over refactoring in this case
    """
    if transmuteObject.moiety.upper() == "NUCLEOSIDE" : 
        ## Read pdb and convert to pdbobject
        nucleosidePdb = UT.TransmuteToJsonNucleoside(transmuteObject.pdb_fname)
        nucleosidePdb.pdb_for_attributes()


        ## Initialise the main dictionary
        molecule = {}

        #-------------------------------- PDB PROPERTIES --------------------------------#
        # Initialise pdb properties dictionary
        pdb_properties = {}

        # Coordinates and Shape
        pdb_properties["Coordinates"] = json.dumps(nucleosidePdb.get_array())
        pdb_properties["Shape"] = json.dumps(nucleosidePdb.get_shape_array())

        # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings (escape characters)
        pdb_properties["Atoms"] = json.dumps(nucleosidePdb.get_atoms())

        # Get the element symbol
        pdb_properties["Symbol"] = json.dumps(nucleosidePdb.get_element_symbol())

        molecule["pdb_properties"] = pdb_properties
        #----------------------------------- IDENTITY -----------------------------------#
        # Initialise the identity list
        identity = []
        # full name
        fullname = nucleosidePdb.get_full_name(transmuteObject.chemistry)
        identity.append(fullname)

        # abbreviated name
#        abbr = nucleosidePdbChemistry
        identity.append(transmuteObject.chemistry)
#
        # molecule chemistry, which is often the same as the residue name in the pdb
        molecule_residuename = nucleosidePdb.get_resname()
        identity.append(molecule_residuename)

        # Nucleobase of the nucleic acid
        nucleobaseName = nucleosidePdb.get_nucleobase(transmuteObject.nucleobase)
        identity.append(nucleobaseName)

        # Check if atoms to build by, stated in the TABLES, are present in the pdb
        nucleosidePdb.validate_atomnames_for_building(chemistry=transmuteObject.chemistry,
                                                    nucleobase=transmuteObject.nucleobase
                                                    )
        molecule["identity"] = json.dumps(identity)

        #------------------------------- TORSIONS AND ANGLES -----------------------------#
        # Initialise the dictionary for the dihedrals and the bond angles
        angles = {}

        # Get Bond angles
        bondangles = nucleosidePdb.get_angles(transmuteObject.angles)

        # Get dihedrals
        torsions = nucleosidePdb.get_angles(transmuteObject.dihedrals)

        angles["bond_angles"] = json.dumps(bondangles)
        angles["dihedrals"] = json.dumps(torsions)

        #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
        molecule["angles"] = angles

        #----------------------------- WRITE OUT A JSON FILE ----------------------------#
        # The json dump() method always requires us to dump it to a file in the current directory
        fname = nucleosidePdb.get_output_name(transmuteObject.chemistry, transmuteObject.conformation, transmuteObject.nucleobase)

        # Get Ducque home
        DUCQUEHOME = SD.return_DUCQUEHOME()

        # Write the inputfile to the Ducque home directory
        with open(DUCQUEHOME + "json/" + fname + ".json", "w") as filejson:
            json.dump(molecule, filejson, indent=4)

        SD.print_writing(f"{DUCQUEHOME}json/{fname}.json")


    elif transmuteObject.moiety.upper() == "LINKER" : 
        ## Read pdb and convert to pdbobject
        linkerPdb = UT.TransmuteToJsonLinker(transmuteObject.pdb_fname)
        linkerPdb.pdb_for_attributes()


        ## Initialise the main dictionary
        molecule = {}

        #-------------------------------- PDB PROPERTIES --------------------------------#
        # Initialise pdb properties dictionary
        pdb_properties = {}

        # Coordinates and Shape
        pdb_properties["Coordinates"] = json.dumps(linkerPdb.get_array())
        pdb_properties["Shape"] = json.dumps(linkerPdb.get_shape_array())

        # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings (escape characters)
        pdb_properties["Atoms"] = json.dumps(linkerPdb.get_atoms())

        # Get the element symbol
        pdb_properties["Symbol"] = json.dumps(linkerPdb.get_element_symbol())

        molecule["pdb_properties"] = pdb_properties
        #----------------------------------- IDENTITY -----------------------------------#
        # Initialise the identity list
        identity = []
        # full name
        fullname = linkerPdb.get_full_name(transmuteObject.chemistry)
        identity.append(fullname)

        # abbreviated name
#        abbr = linkerPdbChemistry
        identity.append(transmuteObject.chemistry)

        # Check if atoms to build by, stated in the TABLES, are present in the pdb
        linkerPdb.validate_atomnames_for_building(chemistry=transmuteObject.chemistry)

        molecule["identity"] = json.dumps(identity)

        #------------------------------- TORSIONS AND ANGLES -----------------------------#
        # Initialise the dictionary for the dihedrals and the bond angles
        angles = {}

        # Get Bond angles
        bondangles = linkerPdb.get_angles(transmuteObject.angles)

        # Get dihedrals
        torsions = linkerPdb.get_angles(transmuteObject.dihedrals)

        angles["bond_angles"] = json.dumps(bondangles)
        angles["dihedrals"] = json.dumps(torsions)

        #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
        molecule["angles"] = angles

        #----------------------------- WRITE OUT A JSON FILE ----------------------------#
        # The json dump() method always requires us to dump it to a file in the current directory
        fname = linkerPdb.get_output_name(transmuteObject.chemistry, transmuteObject.conformation)

        # Get Ducque home
        DUCQUEHOME = SD.return_DUCQUEHOME()

        # Write the inputfile to the Ducque home directory
        with open(DUCQUEHOME + "json/" + fname + ".json", "w") as filejson:
            json.dump(molecule, filejson, indent=4)

        SD.print_writing(f"{DUCQUEHOME}json/{fname}.json")




    elif transmuteObject.moiety.upper() == "NUCLEOBASE" : 
        ## Read pdb and convert to pdbobject
        nucleobasePdb = UT.TransmuteToJsonNucleobase(transmuteObject.pdb_fname)
        nucleobasePdb.pdb_for_attributes()


        ## Initialise the main dictionary
        molecule = {}

        #-------------------------------- PDB PROPERTIES --------------------------------#
        # Initialise pdb properties dictionary
        pdb_properties = {}

        # Coordinates and Shape
        pdb_properties["Coordinates"] = json.dumps(nucleobasePdb.get_array())
        pdb_properties["Shape"] = json.dumps(nucleobasePdb.get_shape_array())

        # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings (escape characters)
        pdb_properties["Atoms"] = json.dumps(nucleobasePdb.get_atoms())

        # Get the element symbol
        pdb_properties["Symbol"] = json.dumps(nucleobasePdb.get_element_symbol())

        molecule["pdb_properties"] = pdb_properties
        #----------------------------------- IDENTITY -----------------------------------#
        # Initialise the identity list
        identity = nucleobasePdb.get_full_name(transmuteObject.chemistry)

        # abbreviated name
#        abbr = nucleobasePdbChemistry
#        identity.append(transmuteObject.chemistry)

        # Check if atoms to build by, stated in the TABLES, are present in the pdb
        nucleobasePdb.validate_atomnames_for_building(transmuteObject.atoms)
        molecule["identity"] = json.dumps(identity)

        #----------------------------- ATOMS FOR ROTATION ----------------------------#
        atoms_rotation = nucleobasePdb.get_atoms_for_rotation(transmuteObject.atoms)
        molecule["atoms_rotation"] = json.dumps(atoms_rotation)

        #----------------------------- WRITE OUT A JSON FILE ----------------------------#
        # The json dump() method always requires us to dump it to a file in the current directory
        fname = nucleobasePdb.get_output_name(transmuteObject.chemistry)

        # Get Ducque home
        DUCQUEHOME = SD.return_DUCQUEHOME()

        # Write the inputfile to the Ducque home directory
        with open(DUCQUEHOME + "json/" + fname + ".json", "w") as filejson:
            json.dump(molecule, filejson, indent=4)

        SD.print_writing(f"{DUCQUEHOME}json/{fname}.json")

    else :
        SD.print_invalid_argument(transmuteObject.moiety, "--moiety")

#    ## Read pdb and convert to pdbobject
#    nucleicAcid = UT.TransmuteToJson(transmuteObject.pdb_fname)
#    nucleicAcid.pdb_for_attributes()
#
#
#    ## Initialise the main dictionary
#    molecule = {}
#
#    #-------------------------------- PDB PROPERTIES --------------------------------#
#    # Initialise pdb properties dictionary
#    pdb_properties = {}
#
#    # Coordinates and Shape
#    pdb_properties["Coordinates"] = json.dumps(nucleicAcid.get_array())
#    pdb_properties["Shape"] = json.dumps(nucleicAcid.get_shape_array())
#
#    # Get Atom namelist             NB: json outputs double quotations as \" XX \" for strings (escape characters)
#    pdb_properties["Atoms"] = json.dumps(nucleicAcid.get_atoms())
#
#    # Get the element symbol
#    pdb_properties["Symbol"] = json.dumps(nucleicAcid.get_element_symbol())
#
#    molecule["pdb_properties"] = pdb_properties
#    #----------------------------------- IDENTITY -----------------------------------#
#    # Initialise the identity list
#    identity = []
#
#    if transmuteObject.moiety.upper() == "NUCLEOSIDE":
#        # full name
#        fullname = nucleicAcid.get_full_name(transmuteObject.chemistry, transmuteObject.moiety)
#
#        # abbreviated name
#        identity.append(transmuteObject.chemistry)
#
#        # molecule chemistry, which is often the same as the residue name in the pdb
#        molecule_residuename = nucleicAcid.get_resname()
#        identity.append(molecule_residuename)
#
#        # Nucleobase of the nucleic acid
#        nucleobaseName = nucleicAcid.get_nucleobase(transmuteObject.nucleobase)
#        identity.append(nucleobaseName)
#
#        # Check if atoms to build by, stated in the TABLES, are present in the pdb
#        nucleicAcid.validate_atomnames_for_building(moietyType=transmuteObject.moiety,
#                                                    chemistry=transmuteObject.chemistry,
#                                                    nucleobase=transmuteObject.nucleobase
#                                                    )
#        molecule["identity"] = json.dumps(identity)
#
#    if transmuteObject.moiety.upper() == "LINKER" or transmuteObject.moiety.upper() == "NUCLEOBASE":
#
#        # full name
#        fullname = nucleicAcid.get_full_name(transmuteObject.chemistry, transmuteObject.moiety)
#        identity.append(fullname)
#
#        # abbreviated name
#        abbr = nucleicAcidChemistry
#        identity.append(transmuteObject.chemistry)
#
#        # Check if atoms to build by, stated in the TABLES, are present in the pdb
#        nucleicAcid.validate_atomnames_for_building(moietyType=transmuteObject.moiety, chemistry=transmuteObject.chemistry)
#
#        molecule["identity"] = json.dumps(identity)
#
#
    #------------------------------- TORSIONS AND ANGLES -----------------------------#
#    if transmuteObject.moiety.upper() == "LINKER" or transmuteObject.moiety.upper() == "NUCLEOSIDE":
#        # Initialise the dictionary for the dihedrals and the bond angles
#        angles = {}
#
#        # Get Bond angles
#        bondangles = nucleicAcid.get_angles(transmuteObject.moiety, transmuteObject.angles)
#
#        # Get dihedrals
#        torsions = nucleicAcid.get_angles(transmuteObject.moiety, transmuteObject.dihedrals)
#
#        angles["bond_angles"] = json.dumps(bondangles)
#        angles["dihedrals"] = json.dumps(torsions)
#
#        #----------------- DUMP EVERYTHING INTO THE MOLECULE DICTIONARY -----------------#
#        molecule["angles"] = angles
#
#
#    #----------------------------- WRITE OUT A JSON FILE ----------------------------#
#    # The json dump() method always requires us to dump it to a file in the current directory
#    fname = nucleicAcid.get_output_name(transmuteObject.chemistry, transmuteObject.moiety, transmuteObject.conformation, transmuteObject.nucleobase)
#
#    # Get Ducque home
#    DUCQUEHOME = SD.return_DUCQUEHOME()
#
#    # Write the inputfile to the Ducque home directory
#    with open(DUCQUEHOME + "json/" + fname + ".json", "w") as filejson:
#        json.dump(molecule, filejson, indent=4)
#
#    SD.print_writing(f"{DUCQUEHOME}json/{fname}.json")
#
#
#
#
#
#
#
#
#

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
        
