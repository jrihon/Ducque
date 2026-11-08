from pathlib import Path

import pytest


@pytest.fixture
def construct_dna_inputfile() -> Path:
    return Path("tests/test-data/construct/dna-homoduplex.binp")


@pytest.fixture
def construct_rna_inputfile() -> Path:
    return Path("tests/test-data/construct/rna-homoduplex.binp")


@pytest.fixture
def construct_dna_rna_inputfile() -> Path:
    return Path("tests/test-data/construct/dna-rna-heteroduplex.binp")
