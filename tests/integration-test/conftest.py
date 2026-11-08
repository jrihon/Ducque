from pathlib import Path

import pytest


@pytest.fixture
def snapshot_main_directory() -> Path:
    return Path("tests/test-data/snapshots")


@pytest.fixture
def nucleobase_pdb_path() -> Path:
    return Path("tests/test-data/pdb/adenine.pdb")


@pytest.fixture
def linker_pdb_path() -> Path:
    return Path("tests/test-data/pdb/phosphate_linker.pdb")


@pytest.fixture
def nucleoside_pdb_path() -> Path:
    return Path("tests/test-data/pdb/dna_guanosine_2endo.pdb")


@pytest.fixture
def construct_dna_inputfile() -> Path:
    return Path("tests/test-data/construct/dna-homoduplex.binp").absolute()


@pytest.fixture
def construct_rna_inputfile() -> Path:
    return Path("tests/test-data/construct/rna-homoduplex.binp")


@pytest.fixture
def randomise_inputfile() -> Path:
    return Path("tests/test-data/randomise/randomise-sequence.rinp")


@pytest.fixture
def transmute_nucleoside_inputfile() -> Path:
    return Path("tests/test-data/transmute/transmute-nucleoside.tinp")


@pytest.fixture
def transmute_linker_inputfile() -> Path:
    return Path("tests/test-data/transmute/transmute-linker.tinp")


@pytest.fixture
def transmute_nucleobase_inputfile() -> Path:
    return Path("tests/test-data/transmute/transmute-nucleobase.tinp")
