import random
from typing import Union

import process_CLI_inputs
from builder.parse_or_write import return_chemistrycode
from builder.builder_library import codex_acidum_nucleicum

CODEX = codex_acidum_nucleicum


def join_chemistry_with_sequence(chemistry : str, sequence) -> str:
    """ For the lazy people, concatenate the chemistry with the sequence they want. Only works for single prompted chemistries"""

    sequence = list(map(lambda x: x.strip(","), sequence))

    nucID = process_CLI_inputs.check_if_chemistry_is_valid(chemistry)

    # Get the abbreviated chemistry code
    abbrChemistry = return_chemistrycode(nucID)

    return ", ".join(abbrChemistry + _base for _base in sequence)


def randomise_sequence(chemistry : str, lengthSequence : int) -> str:
    """ depends on the prompted chemistry and the length of the sequence """

    nucID = process_CLI_inputs.check_if_chemistry_is_valid(chemistry)
    print(nucID)
    abbrChemistry = return_chemistrycode(nucID)

    allCodexKeys = list(CODEX.keys())

    listOfPossibleNucleotides = []
    for i, _code in enumerate(allCodexKeys):
        if _code[:-1] == abbrChemistry:
            listOfPossibleNucleotides.append(allCodexKeys[i])

    return ", ".join(random.choice(listOfPossibleNucleotides) for _ in range(lengthSequence))


def randomise_chemistry(chemistry : list, sequence : list) -> str:
    """ randomise the chemistry for a given sequence """

    # Get the given sequence and chemistry into a list without any comma's
    sequence = list(map(lambda x: x.strip(","), sequence))
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    # Include the different chemistries by parsing from the dictionary, taking the first option and then cutting out the base part. leaving only the chemistry
    listOfChemistries = []
    for _chem in chemistry:
        nucID = process_CLI_inputs.check_if_chemistry_is_valid(_chem)
        abbrChemistry = return_chemistrycode(nucID)
        listOfChemistries.append(abbrChemistry)


    # Concatenate strings by adding a random chemistry, from the prompted list, 
    tmpSequence = []
    for nucleotide in sequence:
        _nucl = random.choice(listOfChemistries) + nucleotide
        tmpSequence.append(_nucl)

    return ", ".join(tmpSequence)


def randomise_sequence_and_chemistry(chemistry : list, lengthSequence : int) -> str:
    """ randomise both the sequence and for a set of given chemistry """
    # Get the given chemistry without the comma's in the strings
    chemistry = list(map(lambda x: x.strip(","), chemistry))

    # Get all the abbreviated chemistry codes from the prompted chemistries you want to randomise
    listOfChemistries = []
    for _chem in chemistry:
        nucID = process_CLI_inputs.check_if_chemistry_is_valid(_chem)
        abbrChemistry = return_chemistrycode(nucID)
        listOfChemistries.append(abbrChemistry)

    allCodexKeys = list(CODEX.keys())

    listOfPossibleNucleotides = []
    for i, _code in enumerate(allCodexKeys):
        _code = _code[:-1]

        if _code in listOfChemistries:
            listOfPossibleNucleotides.append(allCodexKeys[i])


    return ", ".join(random.choice(listOfPossibleNucleotides) for _ in range(lengthSequence))


def write_out_complementary_sequence(compl_seq : Union[str, list]) -> str:
    """ If it is a list, parse the list correctly and output it correctly.
        If it is just a string, output it as as it is now. """

    if isinstance(compl_seq, list):
        compl_seq = list(map(lambda x: x.strip(","), compl_seq))
        output_sequence = ", ".join(compl_seq)
        return output_sequence

    return compl_seq
