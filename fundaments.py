import labyrinth

"""When the user prompts the wrong values or flags, this python script intercepts most errors that happen at the start """

class InputExclusivity(Exception):
    Except = " These flags are mutually exclusive; --pdb     --json "



def check_if_nucleotides_are_valid(input_sequence : list) -> bool:
    """ Check if any of the prompted nucleotides is not valid. """
    # Retrieve the keys of the dictionary from which we parse the data
    keys_of_dict = labyrinth.codex_acidum_nucleicum.keys()
    # Check if any of the prompted nucleotides is not found in the sequence
    report = True

    for NA in input_sequence:
        if NA not in keys_of_dict:
            print("One or more of the nucleotides in the given sequence is invalid. Please check your input file : " + NA + "\n\n")
            report = False
            return report

    return report

