import systemsDucque as SD
from ducquelib.library import TABLE_NUCLEOBASE_MODS


#--nucleobase 
#-position 15.A  (required, where chain is optional)
#-mod Psi        (required)
#-resname pU     (optional, if omitted take the resname of original pdb)
#-orientation HG (optional)
class Nucleobase: 

    def __init__(self, nb_instance: int) -> None:

        self.position = ("", "")
        self.from_mod = ""
        self.to_mod = ""
        self.new_resname = ""
        self.reorient = False
        self.first_match = -1 
        self.nb_instance = nb_instance

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


                    for idx, linePdbContent in enumerate(pdbContent) :
                        if linePdbContent.startswith("ATOM"):
                            if linePdbContent[21] == chainLetter and linePdbContent[22:26].strip() == residueNumber : 
                                IsQueryMatched = True
                                self.first_match = idx
                                break

                if not IsQueryMatched : 
                    SD.print_invalid_argument_nb(position, "position", self.nb_instance)

                self.position = (residueNumber, chainLetter)


            # only a residueNumber was queried, in which case we will
            # residueNumber is then just the variable `position`
            else : 
                # See if the position in the pdb matches with the query `-position`
                IsQueryMatched = False
                with open(pdb_nbase_fname, "r") as pdbFile : 
                    pdbContent = pdbFile.readlines()

                    for idx, linePdbContent in enumerate(pdbContent) :
                        if linePdbContent.startswith("ATOM"):
                            if linePdbContent[22:26].strip() == position : 
                                IsQueryMatched = True
                                self.first_match = idx
                                break

                if not IsQueryMatched : 
                    SD.print_invalid_argument_nb(position, "position", self.nb_instance)
            
                self.position = (position, "")



        # If the self.position has already been set
        else : 
            SD.print_already_set("-position", self.position, position)
            SD.exit_Ducque()

    def set_mod(self, mod): 

        if self.from_mod == "" or self.to_mod == "" : 
            # TODO: Check in the nucleobase modification table if this key exists
            # TODO: Make the nucleobase table modification (done?)

            # There are two ways to use the -mod flag for modification
            # 1. You convert from one nucleobase to the other, and possibly switch its orientation
            # 2. You keep the natural or modified nucleobase and just switch its orientation


            ## 1. modify nucleobase
            if "," in mod: #splittable = yes

                m = mod.split(",")
                try : 
                    from_mod = m[0] 
                    to_mod = m[1] 
                except IndexError : 
                    SD.print_empty_query_nb("-mod", self.nb_instance)
                    SD.exit_Ducque()

                if to_mod == "": #if this variable is empty, the query was incomplete
                    SD.print_invalid_argument_nb(mod, "-mod", self.nb_instance)

                # set attributes
                self.from_mod = from_mod    # from modification
                self.to_mod = to_mod        # to modification 
                
            ## 2. State where we don't do anything, but populate the self.to_mod and self.from_mod
            # Very likely this position will get a reorientation, but that is handled explicitly elsewhere 
            else : 
                self.from_mod = mod     # from modification
                self.to_mod = mod       # to modification


            # Check if self.from_mod and self.to_mod are valid keys in the available nucleobase present in the Ducque library
            try : 
                _ = TABLE_NUCLEOBASE_MODS[self.from_mod]
            except KeyError : 
                SD.print_invalid_key(self.from_mod, "TABLE_NUCLEOBASE_MODS")
                SD.exit_Ducque()

            try : 
                _ = TABLE_NUCLEOBASE_MODS[self.to_mod]
            except KeyError : 
                SD.print_invalid_key(self.to_mod, "TABLE_NUCLEOBASE_MODS")
                SD.exit_Ducque()

        else : 
            SD.print_already_set("-mod", self.from_mod + "," + self.to_mod, mod)
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
        At this point in the software, we have already checked if self.position and self.from_mod are filled in

        There are two options here : 
            1. The `-mod` query only passes a parameter to `self.from_mod`
                If this is the case, this means we anticipate for a simple reorientation, and the residue is not altered chemically.
                We query the residue name from the pdb file
                https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
                residue name 18:20 -> index([17:20])

            2. The `-mod` query contains both a value for `self.from_mod` and `self.to_mod`
                If this is the case, this means that the `-resname` is required, because we modify the nucleobase here.
                We need to know the new name for the molecule
        """

        # case where resname can be inferred from the pdb
        if self.from_mod == self.to_mod : 
        
            if self.new_resname == "" : 
                # We have already checked if there is a match or not, so this is safe
                with open(pdb_nbase_fname, "r") as pdbFile : 
                    matchedLine = pdbFile.readlines()[self.first_match]
                    self.new_resname = matchedLine[17:20].strip()
                    print(self.new_resname)

        # case where resname is mandatory
        elif not self.from_mod == "" and not self.to_mod == "": 
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
        (required) from_mod 
        (optional) to_mod           -> sometimes we do not want a modification, just a reorientation
        (optional) new_resname      -> if not set by the user, default to the current one in the pdb
        (optional) reorient         -> flipping the orientation is not always the goal
        """

        if self.position[0] == "" : 
            SD.print_empty_query_nb('-position', self.nb_instance)
            SD.exit_Ducque()

        if self.from_mod == "" : 
            SD.print_empty_query_nb('-mod', self.nb_instance)
            SD.exit_Ducque()

