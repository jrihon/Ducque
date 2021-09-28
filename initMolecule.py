import json
from numpy import pi, asarray

""" Initialise the molecules (nucleoside or linker) as json objects. These objects are used to parse the data of the nucleosides of interest """

class Nucleoside:
    """ nucleoside = Nucleoside(codex_acidum_nucleicum[nucleic_acid-string][0]) """

    def __init__(self, jsonfile):

        with open(jsonfile, "r") as jsonf:
            self.jsonObject = json.load(jsonf)

        self.array = asarray(json.loads(self.jsonObject["pdb_properties"]["Coordinates"]), dtype=float)
        self.atom_list = json.loads(self.jsonObject["pdb_properties"]["Atoms"])
        self.mol_length = int(json.loads(self.jsonObject["pdb_properties"]["Shape"])[0])

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


class Desmos(Nucleoside):
    """  We can just simply pass this in here for now, since we essentially copy the parent class """

    def get_nucleic_acid_code(self) -> str:
        """ returns the type of chemistry of the nucleic acid. """
        return json.loads(self.jsonObject["identity"])[0]


