from pathlib import Path

from ducque.builder.utils_builder import generate_complementary_sequence
from ducque.refactorlater.process_cli_inputs import build


def test_preprocessing_dna_construct(construct_dna_inputfile: Path) -> None:

    leading_strand, complement, output_filename = build(construct_dna_inputfile)

    assert leading_strand == ["DT", "DC", "DA", "DA", "DC", "DG", "DG", "DT", "DA"]
    assert complement == "HOMO"
    assert output_filename == "dna-homoduplex"


def test_preprocessing_rna_construct(construct_rna_inputfile: Path) -> None:

    leading_strand, complementary_strand, output_filename = build(
        construct_rna_inputfile
    )
    assert leading_strand == ["RU", "RC", "RA", "RA", "RC", "RG", "RG", "RU", "RA"]
    assert complementary_strand == "HOMO"
    assert output_filename == "rna-homoduplex"


def test_generate_complement_dna(construct_dna_inputfile: Path) -> None:
    leading_strand, complementary_strand, _ = build(construct_dna_inputfile)

    # Assert the instancing of the complementary strand
    complementary_strand = generate_complementary_sequence(
        leading_strand, complementary_strand
    )

    assert complementary_strand == [
        "DA",
        "DG",
        "DT",
        "DT",
        "DG",
        "DC",
        "DC",
        "DA",
        "DT",
    ]


def test_generate_complement_rna(construct_rna_inputfile: Path) -> None:
    leading_strand, complementary_strand, _ = build(construct_rna_inputfile)

    # Assert the instancing of the complementary strand
    complementary_strand = generate_complementary_sequence(
        leading_strand, complementary_strand
    )

    assert complementary_strand == [
        "RA",
        "RG",
        "RU",
        "RU",
        "RG",
        "RC",
        "RC",
        "RA",
        "RU",
    ]


def test_generate_dna_with_rna_complement(construct_dna_rna_inputfile: Path) -> None:
    leading_strand, complementary_strand, _ = build(construct_dna_rna_inputfile)

    # Assert the instancing of the complementary strand
    complementary_strand = generate_complementary_sequence(
        leading_strand, complementary_strand
    )

    assert complementary_strand == [
        "RA",
        "RG",
        "RU",
        "RU",
        "RG",
        "RC",
        "RC",
        "RA",
        "RU",
    ]
