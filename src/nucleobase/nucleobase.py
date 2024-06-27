# first file of the WIP
from os.path import isfile

import nucleobase.utils_nucleobase as UN
from nucleobase.utils_nucleobase import NucleobaseMod, PdbFragment
import systemsDucque as SD
from ducquelib.library import TABLE_NUCLEOBASE_MODS
from initMolecule import Nucleobase
from builder import mathematics
from builder import parse_or_write



def modify_nucleobases(PDB_NBASE_FNAME: str, LIST_OF_NBASE_MODIFICATIONS : list[UN.NucleobaseMod]) -> None : 
    """ 
    take in the CLI inputs
    -> pdb file 
        
    -> pre-processed CLI inputs into an Object
        --nucleobase 
        -position 15.A  (required, where chain is optional)
        -mod Psi        (required)
        -resname pU     (optional, takes on resname of the current residue in the original pdb)
        -reorient       (optional)

    return : Outputs a new pdb file 
    """

    # Has already been checked, but double check again won't hurt
    if not isfile(PDB_NBASE_FNAME): 
        SD.print_filenotfound(PDB_NBASE_FNAME)
        SD.exit_Ducque()

    # Read in the file once and process is as many times as needed
    with open(PDB_NBASE_FNAME, "r") as pdbFile : 
        pdbFileContent = pdbFile.readlines()


    # Sort LIST_OF_NBASE_MODIFICATIONS by the start_match attribute
    # This way we can iterate over the LIST_OF_NBASE_MODIFICATIONS later on when writing out the pdb
    sortedListOfModifications: list[NucleobaseMod] = list()
    _tmplist = [x.start_match for x in LIST_OF_NBASE_MODIFICATIONS]
    _tmplist.sort(key=int)

    for start_match in _tmplist : 
        for element in LIST_OF_NBASE_MODIFICATIONS : 
            if start_match == element.start_match : 
                sortedListOfModifications.append(element)



    # Make a list of the PdbFragments
    # A PdbFragment is a class that contains all the information of the residue (from the prompted pdb) we want to modify 
    pdbFragments: list[PdbFragment] = list()

    for modification in LIST_OF_NBASE_MODIFICATIONS : 
        
        # Pass fileContent starting within a precomputed range ; one full residue
        # Append the fragments to the list
#        print(pdbFileContent[modification.start_match:modification.end_match])
        pdbFragments.append( 
                            PdbFragment( pdbFileContent[modification.start_match:modification.end_match] ) 
                            ) 


    # Iterate over both objects, as every ith fragment corresponds the ith modification 
    modifiedPdbFragments : list[PdbFragment] = list()
    for fragment, modification in zip(pdbFragments, LIST_OF_NBASE_MODIFICATIONS): 

        # Query from modifications we need to modify from
        nucleobaseModFromJson = TABLE_NUCLEOBASE_MODS[modification.modify_from][0]
        # A Nucleobase() class retrieves and stores all data from the nucleobase.json file
        nucleobaseFrom = Nucleobase(nucleobaseModFromJson)

        # Check if atoms, to be used rotations, are valid
        # the method prints a valid error message
        if not fragment.validate_atomnames_from_query(nucleobaseFrom):
            SD.exit_Ducque()


        # Query the modification we need to modify to
        nucleobaseModToJson = TABLE_NUCLEOBASE_MODS[modification.modify_to][0]
        nucleobaseTo = Nucleobase(nucleobaseModToJson)




        ## GET THE QUATERNION TO ROTATE WITH

        # get indices of the Coordinates in the array
        indicesArrayFrom = fragment.retrieve_atom_index_MULTIPLE(nucleobaseFrom.atomsRotation)
        indicesArrayTo = parse_or_write.retrieve_atom_index_MULTIPLE(nucleobaseTo, nucleobaseTo.atomsRotation)


        # Generate direction axis for both nucleobases
        # Suppose atomsRotation are [N9, C4, C8]
        # Then we need v1(N9 -> C4) and v2(N9 -> C8)

        # Get plane of the original nucleobase, which is the plane we rotate onto
        directionAxisOriginal = mathematics.get_normal_vector_of_plane( fragment.array[indicesArrayFrom[1]] - fragment.array[indicesArrayFrom[0]],
                                                                          fragment.array[indicesArrayFrom[2]] - fragment.array[indicesArrayFrom[0]]
                                                                            )
        directionAxisModification = mathematics.get_normal_vector_of_plane( nucleobaseTo.array[indicesArrayTo[1]] - nucleobaseTo.array[indicesArrayTo[0]],
                                                                            nucleobaseTo.array[indicesArrayTo[2]] - nucleobaseTo.array[indicesArrayTo[0]]
                                                                            )

        quaternionNucleobase = mathematics.get_quaternion(vector_to_rotate_onto=directionAxisOriginal,
                                                          vector_to_rotate_from=directionAxisModification
                                                 )

        ## Perform actual rotation on nucleobase we want to modify with
        # Move to origin first
        movedToOriginModification = mathematics.move_vector_to_origin(nucleobaseTo.array,
                                                                      distance_to_origin=nucleobaseTo.array[indicesArrayTo[0]]
                                                                      )
        # Rotate accordingly
        rotatedModification = mathematics.rotate_with_quaternion(quaternionNucleobase, movedToOriginModification)

        # Second rotation involves aligned the modified nucleobase with the original nucleobase
        quaternionAlign = mathematics.get_quaternion(vector_to_rotate_from= rotatedModification[indicesArrayTo[1]] - rotatedModification[indicesArrayTo[0]],
                                                     vector_to_rotate_onto= fragment.array[indicesArrayFrom[1]] - fragment.array[indicesArrayFrom[0]] 
                                                    )
        finalRotatedModification = mathematics.rotate_with_quaternion(quaternionAlign, rotatedModification)

        # Move to location of original nucleobase we want to modifiy
        # Move array to Nucleobase object
        nucleobaseTo.array = mathematics.move_vector_to_loc(finalRotatedModification, 
                                                            distance_to_loc=fragment.array[indicesArrayFrom[0]]
                                                            )
#        nucleobaseTo.array = mathematics.move_vector_to_loc(rotatedModification, 
#                                                            distance_to_loc=fragment.array[indicesArrayFrom[0]]
#                                                            )


        # If the user wants to reorient (from watson crick franklin -> hoogsteen), then we invert the plane (directionAxis) we rotate onto
        if modification.reorient : 
            print("Will implement later")

        # Modifiy the original residues with the modifications
        modifiedPdbFragment = UN.change_pdbfragments_by_modification(fragment, nucleobaseFrom, nucleobaseTo)

        # Append to the list
        modifiedPdbFragments.append(modifiedPdbFragment)

    # Write out a pdb file with the modified nucleobases
    UN.write_out_pdb_file_with_modified_nucleobases(PDB_NBASE_FNAME, pdbFileContent, LIST_OF_NBASE_MODIFICATIONS, modifiedPdbFragments) 
