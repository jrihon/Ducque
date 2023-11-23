import sys


""" Contains all the dictionaries and the tools to help the Transmutation work properly.    """

# Identify the nucleic acid chemistry
TABLE_CHEMISTRY = {
        "DNA" : "Deoxyribonucleic acid",
        "RNA" : "Ribonucleic acid",
        "B-HOMODNA" : "beta-Homo Deoxyribonucleic acid",
        "HNA" : "Hexitol Nucleic acid",
        "XYLO" : "Xylose Nucleic acid",
        "2-OME-RNA" : "2-O-Methyl Ribonucleic acid",
        "CENA" : "Cylcohexenyl Nucleic acid",
        "DXYLO" : "deoxy Xylose Nucleic Acid",
        "MNA" : "Morpholino Nucleic Acid",
        "TNA" : "Threose Nucleic Acid",
        "PHOTNA" : "Phosphonate Threose Nucleic Acid",
        }

# the nucleoside moiety accounts for a binding between Phosphate - Carbon
TABLE_LINKER = {
        "DNA" : "Phosphate",
        "RNA" : "Phosphate",
        "B-HOMODNA" : "Phosphate",
        "HNA" : "Phosphate",
        "XYLO" : "Phosphate",
        "2-OME-RNA" : "Phosphate",
        "CENA" : "Phosphate",
        "DXYLO" : "Phosphate",
        "MNA" : "Phosphate",
        "TNA" : "Phosphate",
        "PHOTNA" : "Phosphate", 
        }


# Used for the filename and to fill out the identity dict in the json file
TABLE_NUCLEOBASE = {
        "A" : "Adenosine",
        "C" : "Cytidine",
        "G" : "Guanosine",
        "T" : "Thymidine",
        "U" : "Uracil",
        }



def get_base_type(base : str) -> str:
    """ Retrieve the type of base we will calculate with 
        At the moment, this function is unused, but might be used later when the modified nucleobases are implemented
    """

    if base == "A":
        return "purine"

    if base == "T":
        return "purine"

    if base == "C":
        return "pyrimidine"

    if base == "G":
        return "pyrimidine"

    if base == "U":
        return "pyrimidine"

    # If you've reached this part, then something has gone wrong and it is not clear what the base is
    if True:
        print("It is not clear what the base is. Please check the Residue Name column in the prompted pdb file.\n"
                + "NB: The last character of the string should end as the identifier of one of the five canonical bases.")
        sys.exit(0)





#------------------------ HOLDING ON TO THIS IDEA FOR LATER (?) ----------------------------------------------------
# The general backbone, whether it is a purine or a pyrimidine. This is a nested dictionary.
# Used to parse the indices of the vectors in the molecule's array
#dihedral_dict = {
#        "DNA" : {"backbone" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
#                 "purine" : ["O4'", "C1'", "N9", "C4"],
#                 "pyrimidine" : ["O4'", "C1'","N1", "C2"]
#                 },
#        "RNA" : {"backbone" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
#                 "purine" : ["O4'", "C1'", "N9", "C4"],
#                 "pyrimidine" : ["O4'", "C1'","N1", "C2"]
#                 },
#        }
#
#
# The general backbone, whether it is a purine or a pyrimidine. This is a nested dictionary.
# Used to parse the indices of the vectors in the molecule's array
#angle_dict = {
#        "DNA" : {"backbone" : [""],
#                 "purine" : ["C1'", "N9", "C4"],
#                 "pyrimidine" : ["C1'","N1", "C2"],
#                 },
#        "RNA" : {"backbone" : [""],
#                 "purine" : ["C1'", "N9", "C4"],
#                 "pyrimidine" : ["C1'","N1", "C2"],
#                 }
#        }
#
