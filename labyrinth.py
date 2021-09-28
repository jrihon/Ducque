import sys, os, json
import numpy as np

import initMolecule

import labyrinth_func as LabF                       # Import all the functions to do with creating the nucleic acid duplex
import labyrinth_func_andthefunkybunch as LFFB      # Import all smaller functions
import labyrinth_func_tools3                        # Import the nucleic acid dictionaries


""" Create dictionary of the filename and their molecule code """

CODEX_AN = labyrinth_func_tools3.codex_acidum_nucleicum   # Creates an object of the nucleic acid dictionary for the leading strand

CODEX_CONF = labyrinth_func_tools3.conformations_codex     # Creates an object of the nucleic acid dictionary for the complementary strand



def Architecture(nucleic_acid_list, complement):
    """ This function contains all the functions to create the nucleic acid duplex.

    First up, a bottom-up approach is applied when building the nucleic acid single strand.

    A complementary strand is generated by fitting the nucleosides and linker elements, in a bottom-up approach, to the leading strand.

    Finally, a pdb is outputted, with the two chains in present. """

    # Reverse the list order, because we build the nucleoside from the bottom up.
    nucleic_acid_list.reverse()

    index_lead = 0   # Initiate the index counter for the leading strand. 

    # Parse the first nucleotide in the list.
    nucleic_acid = nucleic_acid_list[0]
    ## Start the leading strand

    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json). This also parses the correct linker.
    nucleoside = initMolecule.Nucleoside(CODEX_AN[nucleic_acid][0]) ; index_lead += nucleoside.mol_length
    linker = initMolecule.Desmos(CODEX_AN[nucleic_acid][1]) ; index_lead += linker.mol_length

    # Position the linker moiety on the nucleoside.
    leading_strand = LabF.position_phosphate_linker(nucleoside, nucleoside.array, linker)

    num_nucl = 1 # Initiate a building block counter

    # This for loop positions the subsequent nucleotides of the leading strand.
    for NA in range(1, len(nucleic_acid_list)):

        # Import the next nucleoside and create an object
        nextnuc_acid = nucleic_acid_list[NA]
        #nextnuc = initMolecule.Nucleoside(CODEX_AN[nextnuc_acid][0]) ; index_lead += nextnuc.mol_length
        #nextlink = initMolecule.Desmos(CODEX_AN[nextnuc_acid][1]) ; index_lead += nextlink.mol_length

        # Parse the dictionary for the previous nucleotide, to append the next nucleotide onto.
        previous_nucleic_acid = nucleic_acid_list[NA - 1]
        # Parse the correct nucleoside and linker and create them as objects
        prevnuc, prevlink = initMolecule.Nucleoside(CODEX_AN[previous_nucleic_acid][0]), initMolecule.Desmos(CODEX_AN[previous_nucleic_acid][1])

        # Assert which conformation fits the strand better and take that one and use it build the nucleotide
        conformations = CODEX_CONF[nextnuc_acid]
        nextnuc_asserted_conf = LabF.assert_leading_strand_nucleotide_conformation(conformations, prevnuc, prevlink, leading_strand)
        nextnuc = initMolecule.Nucleoside(nextnuc_asserted_conf) ; index_lead += nextnuc.mol_length
        nextlink = initMolecule.Desmos(CODEX_AN[nextnuc_acid][1]) ; index_lead += nextlink.mol_length

        # Position the following nucleoside in the sequence
        next_nucleoSIDE_positioned = LFFB.position_next_nucleoside(nextnuc, prevnuc, prevlink, leading_strand)

        # Orient the linker better if need be
        # next_nucleoTIDE_positioned = LabF.orient_previous_linker_better()

        # Position the following linker to create a nucleotide array
        if not (NA + 1) == len(nucleic_acid_list):
            # If this is not the last nucleoside in the sequence, build a linker onto it.
            next_nucleoTIDE_positioned = LabF.position_phosphate_linker(nextnuc, next_nucleoSIDE_positioned, nextlink)

            # Create a tuple of the two arrays and stack them. This becomes the leading strand upon which we continue appending nucleotides.
            leading_strand = np.vstack((next_nucleoTIDE_positioned, leading_strand))
        else:
            # Leave the object as a nucleoside, since it is the last one in the sequence. Create a tuple of the two arrays to finalise the leading strand.
            leading_strand = np.vstack((next_nucleoSIDE_positioned, leading_strand))
            # Substract the last linker's molecule size, since it is not included at the end of the leading strand. This index_lead is now set correctly to start the complementary strand build.
            index_lead -= nextlink.mol_length

        num_nucl += 1


    ## Generate the complementary strand
    # a list of complementary nucleotides
    compl_nucleic_acid_list = LabF.generate_complementary_sequence(nucleic_acid_list, complement)
    # print(compl_nucleic_acid_list)

    # Import the first two leading strand nucleosides, but as a list
    leading_nucleosides = [CODEX_AN[nucleic_acid_list[0]], CODEX_AN[nucleic_acid_list[1]]]

    # Both compl1 and compl2 are a list of json filenames. These become instanced object in the function 'assert_starting_bases_..._strand()'
    compl1, compl2 = CODEX_CONF[compl_nucleic_acid_list[0]], CODEX_CONF[compl_nucleic_acid_list[1]]
    compl2_linker = initMolecule.Desmos(CODEX_AN[compl_nucleic_acid_list[1]][1])
    complementary_strand, index_lead = LabF.assert_starting_bases_of_complementary_strand(compl1, compl2, compl2_linker, leading_nucleosides, leading_strand, index_lead)

    index_compl = complementary_strand.shape[0]                 # Increment the index_complementary integer

    # Add subsequent complementary nucleotides
    for cNA in range(2, len(compl_nucleic_acid_list)):
        # Parse the information of the leading strand nucleotide we want to fit the complementary nucleoside against
        lead_nextnuc_acid = nucleic_acid_list[cNA]
        lead_nextnuc, lead_nextlinker = initMolecule.Nucleoside(CODEX_AN[lead_nextnuc_acid][0]), initMolecule.Desmos(CODEX_AN[lead_nextnuc_acid][1])
        index_lead -= (lead_nextnuc.mol_length + lead_nextlinker.mol_length)

        # Parse the name of the files of the different conformations of the complementary nucleic acid
        compl_nextnuc_acid = compl_nucleic_acid_list[cNA]
        conformations, compl_nextlinker = CODEX_CONF[compl_nextnuc_acid], initMolecule.Desmos(CODEX_AN[compl_nextnuc_acid][1])

        # Parse the previous complementary nucleic acid
        prev_compl_nextnuc_acid = compl_nucleic_acid_list[cNA - 1]
        prev_compl_nuc, prev_compl_linker = initMolecule.Nucleoside(CODEX_AN[prev_compl_nextnuc_acid][0]), initMolecule.Desmos(CODEX_AN[prev_compl_nextnuc_acid][1])
        index_compl -= (prev_compl_linker.mol_length + prev_compl_nuc.mol_length)

        compl_nextnuc_arr = LabF.assert_possible_base_conformations_and_fit(lead_nextnuc, leading_strand, conformations, compl_nextlinker, complementary_strand,
                                                                                                  prev_compl_nuc, prev_compl_linker, index_lead, index_compl)
        complementary_strand = np.vstack((complementary_strand, compl_nextnuc_arr))

        index_compl += (prev_compl_linker.mol_length + prev_compl_nuc.mol_length) + compl_nextnuc_arr.shape[0]

    #------------------------ CREATE THE PDB THAT GOES WITH ARRAY INPUTTED -------------------#
    LabF.create_PDB_from_array_final(leading_strand, nucleic_acid_list, complementary_strand, compl_nucleic_acid_list)
    print("\nNumber of nucleotides in the duplex :" , num_nucl, "\n")

