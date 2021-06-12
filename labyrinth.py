import sys, os, json
import numpy as np

import labyrinth_func as LabF
import labyrinth_func_tools1 as LFT1
import labyrinth_func_tools2 as LFT2

""" Create dictionary of the filename and their molecule code """

codex_acidum_nucleicum = {
'dT': ['json/dna_thymidine.json', 'json/phosphate_linker.json'],
'dA': ['json/dna_adenosine.json', 'json/phosphate_linker.json'],
'dC': ['json/dna_cytosine.json', 'json/phosphate_linker.json'],
'dG': ['json/dna_guanosine.json', 'json/phosphate_linker.json'],


}

def printexit(*arg):
    """Function to print out any numbers of parameters and then exit the software
    This is for debugging purposes	"""
    print(arg)

    sys.exit(0)


def Architecture(nucleic_acid_list):

    # Reverse the list order, because we build the nucleoside from the bottom up
    nucleic_acid_list.reverse()

    ## Parse the first nucleotide in the list
    nucleic_acid = nucleic_acid_list[0]

    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
    nucleoside = LabF.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

    # Parse dictionary for the correct linker segment
    linker = LabF.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

    # LINKER IS POSITIONED
    leading_strand = LabF.position_phosphate_linker(nucleoside, linker)

    # POSITION THE NEXT NUCLEOTIDE
    # import the next nucleoside
    nextnuc_acid = nucleic_acid_list[1]
    nextnuc = LabF.Nucleoside(codex_acidum_nucleicum[nextnuc_acid][0])

    # Parse the dictionary for the previous nucleotide, to append the next nucleotide onto.
    # Parse the correct nucleoside and linker and create them as objects
    previous_nucleic_acid = nucleic_acid_list[1-1]
    prevnuc, prevlink = LabF.Nucleoside(codex_acidum_nucleicum[previous_nucleic_acid][0]), LabF.Desmos(codex_acidum_nucleicum[previous_nucleic_acid][1])

    # Position the following nucleoside in the sequence
    nextnuc_positioned = LabF.position_next_nucleotide(nextnuc, prevnuc, prevlink, leading_strand)

    # create a tuple of the two arrays and stack them. This becomes the leading strand upon which we continue appending nucleotides.
    leading_strand = np.vstack((nextnuc_positioned, leading_strand))

    ######################### CREATE THE PDB THAT GOES WITH ARRAY INPUTTED ####################
    LabF.create_PDB_from_matrix(leading_strand, nucleic_acid_list)

