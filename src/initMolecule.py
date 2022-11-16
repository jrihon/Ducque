import json
from typing import Union, List
from numpy import pi, asarray, arange, ndarray


from builder.utils_labyrinth import generate_complementary_sequence
import builder.parse_or_write as PARSE

""" Initialise the molecules (nucleoside or linker) as json objects. These objects are used to parse the data of the nucleosides of interest """

class Nucleoside:
    """ nucleoside = Nucleoside(codex_acidum_nucleicum[nucleic_acid-string][0]) """

    def __init__(self, jsonfile):

        with open(jsonfile, "r") as jsonf:
            self.jsonObject = json.load(jsonf)

        self.array = asarray(json.loads(self.jsonObject["pdb_properties"]["Coordinates"]), dtype=float)
        self.atom_list = json.loads(self.jsonObject["pdb_properties"]["Atoms"])
        self.mol_length = int(json.loads(self.jsonObject["pdb_properties"]["Shape"])[0])
        self.filename = jsonfile

    def get_dihedral(self, dihedral : str) -> float:
        """ return dihedral value of the queried dihedral """
        # because the dihedral is still inside a dictionary, we need to load the string (json.loads)
        return float(json.loads(self.jsonObject["angles"]["dihedrals"])[dihedral])

    def get_angle(self, angle : str) -> float:
        """ return angle value of the queried angle. Needs to be converted to radians """
        # because the angle is still inside a dictionary, we need to load the string (json.loads)
        return float(json.loads(self.jsonObject["angles"]["bond_angles"])[angle]) * (pi/180)

    def get_base_denominator(self) -> str:
        """ returns the type of base of the nucleic acid. So if the base is Guanosine, return 'G'. """
        return json.loads(self.jsonObject["identity"])[2][-1]

    def get_nucleic_acid_code(self) -> str:
        """ returns the type of chemistry of the nucleic acid. """
        return json.loads(self.jsonObject["identity"])[1]


class Linker(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_nucleic_acid_code(self) -> str:
        """ returns the type of chemistry of the nucleic acid. """
        return json.loads(self.jsonObject["identity"])[0]





class NucleotideSequence():
    """ Create an object to hold all the information about the nucleic acid sequences being prompted by the user """

    def __init__(self, leading_strand_sequence : List[str], complement : Union[List[str], str]):
        """ Standard __init__ function"""

        self.leadingStrandSequence = leading_strand_sequence[::-1]
        self.complementaryStrandSequence = generate_complementary_sequence(self.leadingStrandSequence, complement)




class PdbInstance():
    """ A class created to remove the slow Pandas library from the Daedalus software
        Daedalus used to have the Pandas library because it made it easy to code writing out a *.pdb file.

        I want to get rid of it to improve the speed of Daedalus itself. Instantiating a dataframe and adding to the columns significantly increases the uptime.

        Reads the name of the file and converts the entire file into a workable dataframe.
        By convention of the columns, below, the dataframe is set up like this.
        __________________________________________________________________________________________________________
            # Clarification of the column headers
        Record name :     line 1 - 6
        Atom number :     line 7 - 11
        Atom name :       line 13 - 16
        alternat. loc.:   line 17      (Avaline and Bvaline, delete Bvaline and then clear the column)
        Residue name:     line 18 - 20
        Chain :           line 22
        Sequence number:  line 23 - 26
        X_coord:          line 31 - 38
        Y_coord:          line 39 - 46
        Z_coord:          line 47 - 54
        Occupancy:        line 55 - 60 (no idea what this is; keep)
        Temp. factor:     line 61 - 66 (no idea what this is either; keep)
        Segment id:       line 73 - 76 (no idea what this is, keep it and clear column)
        Element symbol:   line 77 - 78
        Charge:           line 79 - 80
                                                                                            """

    def __init__(self):
#        self.fnamePDB = fnamePDB        # file name of the pdb to be
        self.RecordName = "ATOM"         # this is just written out every single line
        self.AtomNumber = list()
        self.AtomName = list()
        self.ResidueName = list()
        self.Chain = ""
        self.SequenceNumber = list()
        self.x = list()
        self.y = list()
        self.z = list()
        self.Occupancy = "1.00"         # this is just written out every single line
        self.TempFactor = "0.00"        # this is just written out every single line
        self.ElementSymbol = list()


    def SetAtomNumber(self, nucleotideArray : ndarray, addLength : int = 0)  :
        """ Sets the AtomNumber column.
            If we fill in the complementary strand's AtomNumber column, then the variable startNum will not be zero, but equal to the
            length of the startNum of the leading strand's AtomNumber column. """

        rangeAtoms = arange(start=1 + addLength, stop=nucleotideArray.shape[0] + 1 + addLength, dtype=int)

        self.AtomNumber = rangeAtoms
#        self.AtomNumber = [str(x) for x in rangeAtoms]


    def SetAtomName(self, terminal_atomnames : list, nucleotideSequence : list, strandType : str) :
        """  """

        if strandType == "lead" :
            self.AtomName = terminal_atomnames[0] + PARSE.return_PDB_AtomNames_or_ElementSymbol(nucleotideSequence, "Atoms") + terminal_atomnames[1]
        elif strandType == "complementary" :
            self.AtomName = terminal_atomnames[2] + PARSE.return_PDB_AtomNames_or_ElementSymbol(nucleotideSequence, "Atoms") + terminal_atomnames[3]


    def SetResidueName(self, list_of_leading_sequence : list) :
        """  """
        self.ResidueName = PARSE.return_PDB_Residuename(list_of_leading_sequence)


    def SetChainLetter(self, letter : str):
        """  """
        self.Chain = letter


    def SetSequenceNumber(self, nucleotide_sequence : list, strandType : str) :
        """  """
        if strandType == "lead" :
            self.SequenceNumber = PARSE.return_PDB_Sequence(nucleotide_sequence)
        if strandType == "complementary" :
            self.SequenceNumber = PARSE.return_PDB_Sequence(nucleotide_sequence, len(nucleotide_sequence))



    def SetAtomArray(self, terminalArray : ndarray, nucleotideArray : ndarray, strandType : str) -> ndarray :
        """  """
        if strandType == "lead" :
#            leadingArray = vstack((terminalArray[0], nucleotideArray, terminalArray[1]))
            self.x = list(map(lambda x: "{:.3f}".format(x), nucleotideArray[:,0]))
            self.y = list(map(lambda x: "{:.3f}".format(x), nucleotideArray[:,1]))
            self.z = list(map(lambda x: "{:.3f}".format(x), nucleotideArray[:,2]))

        if strandType == "complementary" :
#            leadingArray = vstack((terminalArray[2], nucleotideArray, terminalArray[3]))
            self.x = list(map(lambda x: "{:.3f}".format(x), nucleotideArray[:,0]))
            self.y = list(map(lambda x: "{:.3f}".format(x), nucleotideArray[:,1]))
            self.z = list(map(lambda x: "{:.3f}".format(x), nucleotideArray[:,2]))


    def SetElementSymbols(self, nucleotideSequence : list) :
        """  """
        self.ElementSymbol = ["H"] + PARSE.return_PDB_AtomNames_or_ElementSymbol(nucleotideSequence, "Symbol") + ["H"]
