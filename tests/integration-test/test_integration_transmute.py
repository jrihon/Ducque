import json
from pathlib import Path

from ducque.refactorlater.process_cli_inputs import transmute
from ducque.transmute.transmute import Transmutation, transmute_to_library
from ducque.transmute.utils_transmute import (
    LinkerTransmute,
    NucleobaseTransmute,
    NucleosideTransmute,
)

# from ducque.randomise.randomise import randomiser
# from ducque.refactorlater.process_cli_inputs import build, randomise


def test_transmute_nucleoside_to_library(
    transmute_nucleoside_inputfile: Path,
    snapshot_main_directory: Path,
    tmp_path,
    snapshot,
) -> None:
    """Test the transmute-nucleoside module."""

    snapshot.snapshot_dir = snapshot_main_directory

    # Let preprocessor create inputs for the main function
    transmutation = transmute(transmute_nucleoside_inputfile)

    # Preprocessing module works
    assert isinstance(transmutation, NucleosideTransmute)

    # Modify filename to end up in pytest tmp directory - file
    test_output_dir = Path(tmp_path)
    output_filename = "nucleoside-transmuted.json"
    output_filename_with_ext = f"{output_filename}.json"
    path_to_output = test_output_dir / output_filename_with_ext

    # Transmute nucleobase
    transmute_to_json, transmute_object, molecule = Transmutation(transmutation)
    transmute_to_library(
        transmute_to_json,
        transmute_object,
        molecule,
        path_to_write_to=test_output_dir,
        filename_to_get=output_filename,
    )

    with open(path_to_output, "r") as test_of:
        test_content = test_of.read()

        # assert snapshot
        snapshot.assert_match(
            test_content,  # output to test with
            output_filename_with_ext,  # filename of snapshot that will be created by pytest-snapshot
        )

    with open(path_to_output, "r") as test_of:
        json_content = json.load(test_of)

        assert json_content["pdb_properties"]
        assert json_content["pdb_properties"]["Coordinates"]
        assert json_content["pdb_properties"]["Shape"]
        assert json_content["pdb_properties"]["Atoms"]
        assert json_content["pdb_properties"]["Symbol"]

        assert json_content["identity"]

        assert json_content["angles"]
        assert json_content["angles"]["bond_angles"]
        assert json_content["angles"]["dihedrals"]


def test_transmute_nucleobase_to_library(
    transmute_nucleobase_inputfile: Path,
    snapshot_main_directory: Path,
    tmp_path,
    snapshot,
) -> None:
    """Test the transmute-nucleoside module."""

    snapshot.snapshot_dir = snapshot_main_directory

    # Let preprocessor create inputs for the main function
    transmutation = transmute(transmute_nucleobase_inputfile)

    # Preprocessing module works
    assert isinstance(transmutation, NucleobaseTransmute)

    # Modify filename to end up in pytest tmp directory - file
    test_output_dir = Path(tmp_path)
    output_filename = "nucleobase-transmuted"
    output_filename_with_ext = f"{output_filename}.json"
    path_to_output = test_output_dir / output_filename_with_ext

    # Transmute nucleobase
    transmute_to_json, transmute_object, molecule = Transmutation(transmutation)
    transmute_to_library(
        transmute_to_json,
        transmute_object,
        molecule,
        path_to_write_to=test_output_dir,
        filename_to_get=output_filename,
    )

    with open(path_to_output, "r") as test_of:
        test_content = test_of.read()

        # assert snapshot
        snapshot.assert_match(
            test_content,  # output to test with
            output_filename_with_ext,  # filename of snapshot that will be created by pytest-snapshot
        )

    with open(path_to_output, "r") as test_of:
        json_content = json.load(test_of)

        assert json_content["pdb_properties"]
        assert json_content["pdb_properties"]["Coordinates"]
        assert json_content["pdb_properties"]["Shape"]
        assert json_content["pdb_properties"]["Atoms"]
        assert json_content["pdb_properties"]["Symbol"]

        assert json_content["identity"]

        assert json_content["atoms_rotation"]


def test_transmute_linker_to_library(
    transmute_linker_inputfile: Path,
    snapshot_main_directory: Path,
    tmp_path,
    snapshot,
) -> None:
    """Test the transmute-nucleoside module."""

    snapshot.snapshot_dir = snapshot_main_directory

    # Let preprocessor create inputs for the main function
    transmutation = transmute(transmute_linker_inputfile)

    # Preprocessing module works
    assert isinstance(transmutation, LinkerTransmute)

    # Modify filename to end up in pytest tmp directory - file
    test_output_dir = Path(tmp_path)
    output_filename = "linker-transmuted"
    output_filename_with_ext = f"{output_filename}.json"
    path_to_output = test_output_dir / output_filename_with_ext

    # Transmute nucleobase
    transmute_to_json, transmute_object, molecule = Transmutation(transmutation)
    transmute_to_library(
        transmute_to_json,
        transmute_object,
        molecule,
        path_to_write_to=test_output_dir,
        filename_to_get=output_filename,
    )

    with open(path_to_output, "r") as test_of:
        test_content = test_of.read()

        # assert snapshot
        snapshot.assert_match(
            test_content,  # output to test with
            output_filename_with_ext,  # filename of snapshot that will be created by pytest-snapshot
        )

    with open(path_to_output, "r") as test_of:
        json_content = json.load(test_of)

        assert json_content["pdb_properties"]
        assert json_content["pdb_properties"]["Coordinates"]
        assert json_content["pdb_properties"]["Shape"]
        assert json_content["pdb_properties"]["Atoms"]
        assert json_content["pdb_properties"]["Symbol"]

        assert json_content["identity"]

        assert json_content["angles"]
        assert json_content["angles"]["bond_angles"]
        assert json_content["angles"]["dihedrals"]
