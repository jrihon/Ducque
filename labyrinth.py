import numpy as np
import json
import sys
import os
import labyrinth_func

""" Create dictionary of the filename and their molecule code """
codex_acidum_nucleicum = {
'DT': ['json/dna_thymidine.json', 'json/phosphate_linker.json']


}


def Architecture(nucleic_acid):

    # Parse dictionary for the correct filename; input(DT) - output(dna_thymidine.json)
    nucleo = labyrinth_func.Nucleoside(codex_acidum_nucleicum[nucleic_acid][0])

    # Parse dictionary for the correct linker segment
    #link = labyrinth_func.Desmos(codex_acidum_nucleicum[nucleic_acid][1])

