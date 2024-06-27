import systemsDucque as SD
from initMolecule import Nucleobase
from ducquelib.library import TABLE_NUCLEOBASE_MODS
from numpy import array, ndarray, asarray, zeros


#--nucleobase 
#-position 15.A  (required, where chain is optional)
#-mod Psi        (required)
#-resname pU     (optional, if omitted take the resname of original pdb)
#-reorient       (optional)
class NucleobaseMod: 

    def __init__(self, nb_instance: int) -> None:
        """ 
            nb_instance: The user queries a variable amount of modifications to apply to the prompted model.
                         The `nb_instance` variable just counts the nth time a modification is queried.

            start_match: stores the index at which the position was first found.
                         This means we start parsing the PdbFragment at position[start_match], 
                            so we do not have to reiterate over and over again
        """

        self.position: tuple[str,str] = ("", "")
        self.modify_from: str = ""
        self.modify_to: str = ""
        self.new_resname: str = ""
        self.reorient: bool = False
        self.start_match: int = -1 
        self.end_match: int = -1 
        self.nb_instance: int = nb_instance

    def set_position(self, position, pdb_nbase_fname): 
        """ 
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        Chain is at position 22 - index(21)
        residueNumber is at position 23-26 - index([22:26])
        """

        if self.position[0] == "" :  

            if "." in position : #splittable = yes
                p = position.split(".")
                try : 
                    residueNumber = p[0]
                    chainLetter = p[1]
                except IndexError : 
                    SD.print_empty_query_nb("-position", self.nb_instance)
                    SD.exit_Ducque()

                if chainLetter == "": #if this variable is empty
                    SD.print_invalid_argument_nb(position, "-position", self.nb_instance)

                # See if the position in the pdb matches with the query `-position`
                IsQueryMatched = False
                with open(pdb_nbase_fname, "r") as pdbFile : 
                    pdbContent = pdbFile.readlines()


                    IsQueryMatched = False
                    for idx, linePdbContent in enumerate(pdbContent) :
                        if linePdbContent.startswith("ATOM"):
                            if linePdbContent[21] == chainLetter and linePdbContent[22:26].strip() == residueNumber : 
                                # If we find the first match, remember the index
                                if not IsQueryMatched :
                                    self.start_match = idx
                                    IsQueryMatched = True
                            else : 
                                # If we match on a residue that does not correspond with our matching residue, we are at a new residue
                                # => remember the idx (exclusive range ending), break the loop as we are finished iterating
                                if IsQueryMatched :
                                    self.end_match = idx 
                                    break

                if IsQueryMatched == 0 : 
                    SD.print_invalid_argument_nb(position, "position", self.nb_instance)

                self.position = (residueNumber, chainLetter)


            # only a residueNumber was queried, in which case we will
            # residueNumber is then just the variable `position`
            else : 
                # See if the position in the pdb matches with the query `-position`
                with open(pdb_nbase_fname, "r") as pdbFile : 
                    pdbContent = pdbFile.readlines()

                    IsQueryMatched = False
                    for idx, linePdbContent in enumerate(pdbContent) :
                        if linePdbContent.startswith("ATOM"):
                            if linePdbContent[22:26].strip() == position : 
                                # If we find the first match, remember the index
                                if not IsQueryMatched :
                                    self.start_match = idx
                                    IsQueryMatched = True
                            else : 
                                # If we match on a residue that does not correspond with our matching residue, we are at a new residue
                                # => remember the idx (exclusive range ending), break the loop as we are finished iterating
                                if IsQueryMatched :
                                    self.end_match = idx 
                                    break

                if not IsQueryMatched : 
                    SD.print_invalid_argument_nb(position, "position", self.nb_instance)
            
                self.position = (position, "")



        # If the self.position has already been set
        else : 
            SD.print_already_set("-position", self.position, position)
            SD.exit_Ducque()

    def set_mod(self, mod): 

        if self.modify_from == "" or self.modify_to == "" : 
            # TODO: Check in the nucleobase modification table if this key exists
            # TODO: Make the nucleobase table modification (done?)

            # There are two ways to use the -mod flag for modification
            # 1. You convert from one nucleobase to the other, and possibly switch its orientation
            # 2. You keep the natural or modified nucleobase and just switch its orientation


            ## 1. modify nucleobase
            if "," in mod: #splittable = yes

                m = mod.split(",")
                try : 
                    modify_from = m[0] 
                    modify_to = m[1] 
                except IndexError : 
                    SD.print_empty_query_nb("-mod", self.nb_instance)
                    SD.exit_Ducque()

                if modify_to == "": #if this variable is empty, the query was incomplete
                    SD.print_invalid_argument_nb(mod, "-mod", self.nb_instance)

                # set attributes
                self.modify_from = modify_from    # from modification
                self.modify_to = modify_to        # to modification 
                
            ## 2. State where we don't do anything, but populate the self.modify_to and self.modify_from
            # Very likely this position will get a reorientation, but that is handled explicitly elsewhere 
            else : 
                self.modify_from = mod     # from modification
                self.modify_to = mod       # to modification


            # Check if self.modify_from and self.modify_to are valid keys in the available nucleobase present in the Ducque library
            try : 
                _ = TABLE_NUCLEOBASE_MODS[self.modify_from]
            except KeyError : 
                SD.print_invalid_key(self.modify_from, "TABLE_NUCLEOBASE_MODS")
                SD.exit_Ducque()

            try : 
                _ = TABLE_NUCLEOBASE_MODS[self.modify_to]
            except KeyError : 
                SD.print_invalid_key(self.modify_to, "TABLE_NUCLEOBASE_MODS")
                SD.exit_Ducque()

        else : 
            SD.print_already_set("-mod", self.modify_from + "," + self.modify_to, mod)
            SD.exit_Ducque()


    def set_new_resname(self, new_resname): 

        # this is up to the user to satisfy their need?
        # or maybe make a custom table for the possible resname options - for the GUI at least

        if self.new_resname == "" : 
            self.new_resname = new_resname
        else : 
            SD.print_already_set("-new_resname", self.new_resname, new_resname)
            SD.exit_Ducque()


    def set_reorient(self): 
        # checking if the parameter is already set or not is going to require more instructions than just overriding a Boolean
        self.reorient = True


    def check_if_resname_queried(self, pdb_nbase_fname: str):
        """
        At this point in the software, we have already checked if self.position and self.modify_from are filled in

        There are two options here : 
            1. The `-mod` query only passes a parameter to `self.modify_from`
                If this is the case, this means we anticipate for a simple reorientation, and the residue is not altered chemically.
                We query the residue name from the pdb file
                https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
                residue name 18:20 -> index([17:20])

            2. The `-mod` query contains both a value for `self.modify_from` and `self.modify_to`
                If this is the case, this means that the `-resname` is required, because we modify the nucleobase here.
                We need to know the new name for the molecule
        """

        # case where resname can be inferred from the pdb
        if self.modify_from == self.modify_to : 
        
            if self.new_resname == "" : 
                # We have already checked if there is a match or not, so this is safe
                with open(pdb_nbase_fname, "r") as pdbFile : 
                    matchedLine = pdbFile.readlines()[self.start_match]
                    self.new_resname = matchedLine[17:20].strip()
                    print(self.new_resname)

        # case where resname is mandatory
        elif not self.modify_from == "" and not self.modify_to == "": 
            if self.new_resname == "" : 
                SD.print_empty_query_nb("-resname", self.nb_instance)
                SD.exit_Ducque()

        else : 
            # something shoddy happened
            SD.exit_Ducque()

        

    def check_empty_attributes(self): 
        """ 
        Check if any of the required attributes are set or not
        (required) position         -> only the first of the tuple is required
        (required) modify_from 
        (optional) modify_to           -> sometimes we do not want a modification, just a reorientation
        (optional) new_resname      -> if not set by the user, default to the current one in the pdb
        (optional) reorient         -> flipping the orientation is not always the goal
        """

        if self.position[0] == "" : 
            SD.print_empty_query_nb('-position', self.nb_instance)
            SD.exit_Ducque()

        if self.modify_from == "" : 
            SD.print_empty_query_nb('-mod', self.nb_instance)
            SD.exit_Ducque()







class PdbFragment: 


    def __init__(self, pdbFragmentContent : list[str], ) -> None:
        """ Receives a list of assigned pdb fragment to convert to a proper object 
            
            What we need from the pdb residue is : 
            - atomname list (list[str])
            - atomic coordinates  (np.ndarray)
            - residue number (int)
            - chain letter (str)
            - residue name (str) -> in case not prompted, we default to the present one

            What we already have from the nucleobase.json file
            - coordinate of moiety
            - shape of moiety 
            - atomname list 
            - element list 
            - name of moiety (identity)
            - atoms to rotate by


            Atom name :       line 13 - 16
            Residue name:     line 18 - 20
            Chain :           line 22
            Sequence number:  line 23 - 26
            X_coord:          line 31 - 38 )
            Y_coord:          line 39 - 46 ) self.array
            Z_coord:          line 47 - 54 )
        """

        # Parse the position from the pdb file (residuenumber, chain)
        self.chainLetter: str = pdbFragmentContent[0][21]
        self.residueName: str = pdbFragmentContent[0][17:20].strip()
        self.residueNumber: int = int(pdbFragmentContent[0][22:26])

        atomNames: list[str] = list()
        xCoords = list()
        yCoords = list()
        zCoords = list()
        elements = list()

        for line in pdbFragmentContent :
#            print(line)
#
#            # Check whenever the residueNumber has been incremented
#            residueNumberUnchanged = self.residueNumber == line[22:26].strip()
#
#            # if incremented, break the loop 
#            if not residueNumberUnchanged : 
#                break

            # parse atomnames
            atomNames.append(line[12:16].strip())

            # parse coordinates
            xCoords.append(line[30:38])
            yCoords.append(line[38:46])
            zCoords.append(line[46:54])

            elements.append(line[76:78].strip())


        self.atomNames : list[str] = atomNames
        self.array: ndarray = array([xCoords, yCoords, zCoords], dtype=float).T
        self.elements = elements



    def validate_atomnames_from_query(self, nucleobase : Nucleobase) -> bool :
        """ 
            modAtomNames is a list of atom names e.g. : [N9, C4, C8]
            The atom names come from the queried pdb file and need to match the `modify_from` queried nucleobase
        """

        for atomRotation in nucleobase.atomsRotation: 
            if atomRotation not in self.atomNames : 
                print(f"[INVALID QUERY]   : Atomname `{atomRotation}`, from `{nucleobase.identity}.json` not found in prompted pdb, at position `{(self.residueNumber, self.chainLetter)}`")
                return False

        return True


    def retrieve_atom_index_MULTIPLE(self, atoms : list) :
        """ Retrieves the index in the json_object.array of the atom of interest
            This integer will be used to retrieve the vector of the atom of interest """
        array_of_indexes = zeros(len(atoms), dtype=int)

        for i in range(len(atoms)):
            array_of_indexes[i] = self.atomNames.index(atoms[i])

        return list(array_of_indexes)


def change_pdbfragments_by_modification(pdbFragment : PdbFragment, nucleobaseFrom : Nucleobase, nucleobaseTo : Nucleobase) -> PdbFragment : 
    """
    From the PdbFragment() class, we need only to modify the 
        `.atomnames` 
        `.array` 
        attributes

    From the Nucleobase() class, we need to read from the 
        `.atom_list`
        `.array`
        `.atomsRotation` 
        attributes

    """

    # We want to retrieve the index of the atoms in the residue we want to modifiy
    # Take the `From` molecule and parse the original PdbFragment
    idx_list: list[int] = list()
#    print(nucleobaseFrom.atom_list)
#    print(pdbFragment.atomNames)
    for atomFrom in nucleobaseFrom.atom_list : 
        for idx, atomPdb in enumerate(pdbFragment.atomNames) :
            if atomFrom == atomPdb : 
                idx_list.append(idx)

    # All atoms from the queried nucleobaseFrom have to be found in the original residue 
    # so we can remove them correctly and insert the modified nucleobase (nucleobaseTo)
    # If the `idx_list` does not match the length of the `pdbFragment.atomNames`, then something went wrong

    if len(idx_list) != len(nucleobaseFrom.atom_list) : 
        SD.print_mismatch_nbase(nucleobaseFrom.identity, pdbFragment)
        SD.exit_Ducque()

    # make new instances of the attributes we need to change
    coordinatesPdb = list()
    atomNamesPdb = list()
    elementsPdb = list()

    # On the first instance of a match, we export the entire nucleobaseTo (the modified nucleobase) to the pdbFragment
    isAlreadyMatched = False
    for ith, (coordinate, atomname, element) in enumerate(zip(pdbFragment.array, pdbFragment.atomNames, pdbFragment.elements)): 

        if ith in idx_list : 
            
            # Skip the original nucleobase atoms
            if isAlreadyMatched : 
                continue
            else : 
                # extend the growing list with all this data
                coordinatesPdb.extend(nucleobaseTo.array)
                atomNamesPdb.extend(nucleobaseTo.atom_list)
                elementsPdb.extend(nucleobaseTo.elements)
                # set to True so we do this just once
                isAlreadyMatched = True

        else : 
            coordinatesPdb.append(coordinate)
            atomNamesPdb.append(atomname)
            elementsPdb.append(element)

    pdbFragment.array = asarray(coordinatesPdb, dtype=float)
    pdbFragment.atomNames = atomNamesPdb
    pdbFragment.elements = elementsPdb

    # This returns a modified version of the pdb fragment
    return pdbFragment




def write_out_pdb_file_with_modified_nucleobases(PDB_NBASE_FNAME: str, pdbFileContent : list[str],
                                                 LIST_OF_NBASE_MODIFICATIONS : list[NucleobaseMod],
                                                 modifiedPdbFragments: list[PdbFragment]
                                                 ) -> None : 
    """ 
    Write out the pdb file with modifications implemented
    """


#    def convert_integer_to_string(number: int) -> str : 
#
#        if number > 10 : 
#            return "    " + str(number)
#        elif number > 100 : 
#            return "  " + str(number)
#        elif number > 1000 : 
#            return " " + str(number)
#        else :  # above 1000
#            return str(number)


    # Change name of the output file
    outputfile = PDB_NBASE_FNAME.split(".")[:-1] # grab everything except the file extension
    outputfile.append("_modified.pdb")
    outputfile = "".join(outputfile)

    # Keep track of the place where our modified-residue should go
    counter = 0
    start_match = LIST_OF_NBASE_MODIFICATIONS[counter].start_match
    end_match = LIST_OF_NBASE_MODIFICATIONS[counter].end_match
    matchFound = False

    # make a new pdb file
    with open(outputfile, "w") as pdbFileOutput :

        for i, line in enumerate(pdbFileContent): 
#            atomNumber = i + 1

            # If an atom-line is encountered
            if line.startswith("ATOM"): 

#                # Make the prefix : 
#                prefix = "ATOM  " + convert_integer_to_string(atomNumber)

                # Check if matchFound has ended; we are passed iterating over the residue we intended to modify 
                if end_match == i :
                    # Set values
                    counter += 1
                    matchFound = False

                    # Possible we index outside of the array
                    try : 
                        start_match = LIST_OF_NBASE_MODIFICATIONS[counter].start_match
                        end_match = LIST_OF_NBASE_MODIFICATIONS[counter].end_match
                    except IndexError: 
                        # Set to zero, because have gone through all modifications and 
                        # the rest is going to get written from the original file
                        start_match = 0
                        end_match = 0


                    pdbFileOutput.write(line) 
                    continue


                # If a match has been found, continue the loop,
                # as the entire residue has been modified
                if matchFound : 
                    continue

                if start_match == i :
                    matchFound = True
                    modifiedFragment = modifiedPdbFragments[counter] # We will never index out-of-bonds here, so no need for (try except)
                    newResidueName = LIST_OF_NBASE_MODIFICATIONS[counter].new_resname

                    # Pump and dump all the data into the pdb file being written
                    for j in range(len(modifiedFragment.atomNames)):
                        line = ["ATOM",
                                (i + j + 1),
                                modifiedFragment.atomNames[j],
                                newResidueName,
                                modifiedFragment.chainLetter,
                                modifiedFragment.residueNumber,
                                # x coordinate              , y coordinate                , z coordinate
                                modifiedFragment.array[j][0], modifiedFragment.array[j][1],modifiedFragment.array[j][2],
                                "1.00", "0.00",
                                modifiedFragment.elements[j] ]

#                        pdbFileOutput.write("%-4s  %5s %4s %3s %s%4d    %8s%8s%8s%6s%6s          %2s\n" % tuple(line))
                        pdbFileOutput.write("%-4s  %5s %4s %3s %s%4d    %8.3f%8.3f%8.3f%6s%6s          %2s\n" % tuple(line))


                # If nothing matches or the curent residue is not one we want to modify, just stream the original content into the new pdb
                else : 
                    pdbFileOutput.write(line) 
            #
            # If a TER-line is encountered
            elif line.startswith("TER"): 
                pdbFileOutput.write(line) 
            #
            # Skip any other lines
            else : 
                continue

#    ATOM      4  H5'  DT J   1       4.851   4.813  33.865  1.00  0.00           H
#    ATOM      5 H5''  DT J   1       6.303   5.832  33.718  1.00  0.00           H
#    ATOM      6  C4'  DT J   1       6.305   4.096  32.467  1.00  0.00           C
#    ATOM      7  H4'  DT J   1       7.036   3.609  33.112  1.00  0.00           H

    # 1. KEYWORD  : ATOM
    # 2. The atom number needs to be actively modified for the entire file.
    # 3. AtomNames come from NucleobaseMod
    # 4. ResidueName is to be altered to the modification (NucleobaseMod.new_resname)
    # 5. ChainLetter stays the same (NucleobaseMod.position[1])
    # 6. ResidueNumber stays the same (NucleobaseMod.position[0])
    # 7. Coordinates are changed at the site itself
    # 8. All the things after are whatever
    # 9. ElementSymbols follow the same as AtomNames (NucleobaseMod)

    ## => 
    # Whenever we output, we can pass the unmodified residues as is, except for their atom number

    # Upon encountering a residue we have modified, we need to do two things : 
    #   Take out the coordinates and atomnames of the original nucleobase
    #   -> All atomnames that do not contain an apostrophe are per definition nucleobase atoms
    #       Instead of stacking a new array, just do not include old coordinates in the new file, 
    #       just include the new ones and skip the old ones 'idx, line in enumerate(xxx): if idx == oldcoordinate: continue'
    # 
    #   Insert the coordinates and atomnames of the modified nucleobase
    #   -> Insert the modified array at the instance of the base-atom (i.e. N9 for purines, N1 for pyrimdines)
    #   -> Insert the new set of atomnames also at that instance 
