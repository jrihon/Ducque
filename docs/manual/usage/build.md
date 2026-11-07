# Build

The **Build** module allows to build a structure with a simple query.
The build inputfiles are highlighted with the `.binp` file extension (**b**uild-**inp**ut).

```bash
Ducque --build INPUTFILE.binp
```

To allow for high customisability of model duplexes, the user can build the following way:

- Leading and complementary strand with the same chemistry (homoduplex) (example A).
- A difference in chemistry between both strands (heteroduplex) (example B).
- Complementary strand with various chemistries (custom antisense strands) (example C).

Examples of possibilities of queries accepted by Ducque:

**A.** Homoduplex

```text title="dna_homo_duplex.binp"
--sequence DT, DC, DA, DA, DC, DG, DG, DT, DA
--complement HOMO
--pdbname dna_homo_duplex
```

**B.** Heteroduplex

```text title="rna_mna_heteroduplex.binp"
--sequence RA, RG, RU, RU, RG, RC, RC, RA, RU
--complement MNA
--pdbname rna_mna_heteroduplex
```

**C.** Custom antisense strand

```text title="rna_2ome-rna_heteroduplex.binp"
--sequence   RC, RG,  RA,  RA, RG, RC, RC,  RA
--complement RG, RC, 2MU, 2MU, RC, RG, RG, 2MU
--pdbname rna_2ome-rna_heteroduplex
```

## Sequence

`--sequence`: A list of comma-separated values. The values are specific nucleic acids that are implemented into Ducque.
A list of values readily available to the user is presented in [Available input by default](available-inputs.md).

## Complement

`--complement`: This can be one of three different types of queries.

1. The `HOMO` input, to denote that the complementary strand should be built with the same chemistry as the leading strand.
   The complementary sequence is inferred from the `--sequence` query.
2. The specific chemistry for the complementary sequence, such as `DNA`, `RNA`, `MNA` etc.
   The complementary sequence is inferred from the `--sequence` query.
3. A list of comma-separated values, as shown in [`--sequence`](#sequence).
   A list of values readily available to the user is presented in [Available input by default](available-inputs.md).
   This feature is only available from the CLI, as this requires manual entry of the values.

## PDB name

`--pdbname`: The name of the **output** pdb-file. This can be appended by the `.pdb` extension, but is not required.
