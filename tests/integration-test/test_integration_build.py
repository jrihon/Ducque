from pathlib import Path

import pytest

from ducque.builder.builder import Architecture
from ducque.refactorlater.process_cli_inputs import build


@pytest.mark.parametrize(
    "construct_inputfile_fixture",
    ["construct_dna_inputfile", "construct_rna_inputfile"],
)
def test_construct_build(
    snapshot_main_directory: Path,
    construct_inputfile_fixture,
    request: pytest.FixtureRequest,
    tmp_path,
    snapshot,
) -> None:
    """Test the construct module."""
    construct_inputfile = request.getfixturevalue(construct_inputfile_fixture)
    snapshot.snapshot_dir = snapshot_main_directory

    # Let preprocessor create inputs for the main function
    leading_strand, complementary_strand, output_filename = build(construct_inputfile)
    output_filename = output_filename + ".pdb"

    # Modify filename to end up in pytest tmp directory - file
    test_outputfile = Path(tmp_path) / output_filename

    # Build model
    Architecture(leading_strand, complementary_strand, str(test_outputfile))

    with open(test_outputfile, "r") as test_of:
        # read content of file
        test_content = test_of.read()

        # assert snapshot
        snapshot.assert_match(
            test_content,  # output to test with
            output_filename,  # filename of snapshot that will be created by pytest-snapshot
        )
