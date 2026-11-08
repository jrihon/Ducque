from pathlib import Path

import pytest

from ducque.randomise.randomise import randomiser
from ducque.refactorlater.process_cli_inputs import build, randomise


@pytest.mark.parametrize(
    "randomise_inputfile_fixture",
    [
        "randomise_homo_inputfile",
        "randomise_hetero_inputfile",
    ],
)
def test_randomiser(
    randomise_inputfile_fixture,
    request: pytest.FixtureRequest,
    tmp_path,
) -> None:
    """Test the randomise module."""
    construct_inputfile = request.getfixturevalue(randomise_inputfile_fixture)

    # Let preprocessor create inputs for the main function
    chemistry, length_sequence, sequence, complement, output_filename = randomise(
        construct_inputfile
    )
    output_filename = output_filename + ".binp"

    # Modify filename to end up in pytest tmp directory - file
    test_outputfile = Path(tmp_path) / output_filename

    # Build randomised sequence
    randomiser(chemistry, length_sequence, sequence, complement, str(test_outputfile))

    # Satisfy three commandline for the output
    with open(test_outputfile, "a") as test_of:
        test_of.write("--pdbname foobar")

    leading_strand, complementary_strand, _ = build(test_outputfile)

    assert len(leading_strand) == length_sequence
    assert complementary_strand == "HOMO"
