import systemsDucque as SD



#--nucleobase 
#-position 15.A  (required, where chain is optional)
#-mod Psi        (required)
#-resname pU     (optional, if omitted take the resname of original pdb)
#-orientation HG (optional)
class Nucleobase: 

    def __init__(self) -> None:

        self.position = ("", "")
        self.from_mod = ""
        self.to_mod = ""
        self.new_resname = ""
        self.reorient = False

    def set_position(self, position, pdb_nbase_fname): 

        if self.position[0] == "" :  

            if "." in position : #splittable = yes
                p = position.split(".")
                try : 
                    residueNumber = p[0]
                    chainLetter = p[1]
                except IndexError : 
                    SD.print_empty_query("-position")
                    SD.exit_Ducque()

                if chainLetter == "": #if this variable is empty
                    SD.print_invalid_argument(position, "-position")

                self.position = (residueNumber, chainLetter)

                # TODO: check if position queried is a valid position against the queried pdb file

            # only a residueNumber was queried, in which case we will
            else : 
                self.position = (position, "")

                # TODO: check if position queried is a valid position against the queried pdb file




        # If the self.position has already been set
        else : 
            SD.print_already_set("-position", self.position, position)
            SD.exit_Ducque()

    def set_mod(self, mod): 

        if self.from_mod == "" or self.to_mod == "" : 
            # TODO: Check in the nucleobase modification table if this key exists
            # TODO: Make the nucleobase table modification (done?)

            # There are two ways to use the --nbase flag for modification
            # 1. You convert from one nucleobase to the other, and possibly switch its orientation
            # 2. You keep the natural or modified nucleobase and just switch its orientation


            # 1. modify nucleobase
            if "," in mod: #splittable = yes

                m = mod.split(",")
                try : 
                    from_mod = m[0]
                    to_mod = m[1]
                except IndexError : 
                    SD.print_empty_query("-mod")
                    SD.exit_Ducque()

                if to_mod == "": #if this variable is empty
                    SD.print_invalid_argument(mod, "-mod")
                
            # 2. State where we don't do anything, but populate the self.to_mod and self.from_mod
            # Very likely this position will get a reorientation, but that is handled explicitly elsewhere 
            else : 

                # from modification 
                self.from_mod = mod

                # to modification 
                self.to_mod = mod

            # Check if self.from_mod and self.to_mod are valid keys in the available nucleobase present in the Ducque library
#
#
#
#
#
#
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

