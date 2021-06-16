import sys, os, json
import numpy as np

#from scipy.spatial.transform import Rotation as R
import labyrinth_func as LabF

""" Create dictionary of the filename and their molecule code """

codex_acidum_nucleicum = {
'dT': ['json/dna_thymidine.json', 'json/phosphate_linker.json'],
'dA': ['json/dna_adenosine.json', 'json/phosphate_linker.json'],
'dC': ['json/dna_cytidine.json', 'json/phosphate_linker.json'],
'dG': ['json/dna_guanosine.json', 'json/phosphate_linker.json'],
}

def Architecture(nucleic_acid_list):

    # Reverse the list order, because we build the nucleoside from the bottom up
    nucleic_acid_list.reverse()

    ## Parse the first nucleotide in the list
    nucleic_acid = nucleic_acid_list[0]

    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
    nucleoside = LabF.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])
#    test_angle = 0.094 * np.pi
#    nucleoside.array = LFT1.rotate_with_quaternion(R.from_quat([0, 0, np.sin(test_angle/2), np.cos(test_angle/2)]), nucleoside.array)
    # Parse dictionary for the correct linker segment
    linker = LabF.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

    # LINKER IS POSITIONED
    leading_strand = LabF.position_phosphate_linker(nucleoside, nucleoside.array, linker)

    for NA in range(1, len(nucleic_acid_list)):
        print("-----------------------------------------------------------")
        # POSITION THE NEXT NUCLEOTIDE
        # import the next nucleoside and create an object
        nextnuc_acid = nucleic_acid_list[NA]
        nextnuc = LabF.Nucleoside(codex_acidum_nucleicum[nextnuc_acid][0])
        nextlink = LabF.Desmos(codex_acidum_nucleicum[nextnuc_acid][1])

#        nextnuc.array = LFT1.rotate_with_quaternion(R.from_quat([0, 0, np.sin(test_angle/2), np.cos(test_angle/2)]), nucleoside.array)

        # Parse the dictionary for the previous nucleotide, to append the next nucleotide onto.
        # Parse the correct nucleoside and linker and create them as objects
        previous_nucleic_acid = nucleic_acid_list[NA - 1]
        prevnuc, prevlink = LabF.Nucleoside(codex_acidum_nucleicum[previous_nucleic_acid][0]), LabF.Desmos(codex_acidum_nucleicum[previous_nucleic_acid][1])

        # Position the following nucleoside in the sequence
        next_nucleoSIDE_positioned = LabF.position_next_nucleoside(nextnuc, prevnuc, prevlink, leading_strand)

        # Position the corresponding linker onto the next_nucleoSIDE
        if not (NA + 1) == len(nucleic_acid_list):
            # Position the following linker to create a nucleotide array
            next_nucleoTIDE_positioned = LabF.position_phosphate_linker(nextnuc, next_nucleoSIDE_positioned, nextlink)

            # Create a tuple of the two arrays and stack them. This becomes the leading strand upon which we continue appending nucleotides.
            leading_strand = np.vstack((next_nucleoTIDE_positioned, leading_strand))
        else:
            # Leave the object as a nucleoside, since it is the last one in the sequence
            # Create a tuple of the two arrays and stack them. This becomes the leading strand upon which we continue appending nucleotides.
            leading_strand = np.vstack((next_nucleoSIDE_positioned, leading_strand))

    ######################### CREATE THE PDB THAT GOES WITH ARRAY INPUTTED ####################
    LabF.create_PDB_from_matrix(leading_strand, nucleic_acid_list)

