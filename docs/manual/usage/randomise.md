# Randomise

This module is implemented to facilitate writing and producing sequences.
The **randomisation** inputfiles are highlighted with the `.rinp` file extension (**r**andomise-**inp**ut).

```bash
Ducque --randomise INPUTFILE.rinp
```

- A variable sequence with a specified chemistry (example A).
- A way to prepend a given sequence with a prompted chemistry.
  This method is not a randomisation, but allows the user to write the desired sequence faster.
  This functionality was more fit to be introduced in this module than anywhere else (example B).

This query outputs a `binp` file that can be readily passed to the [build module](build.md).

Examples of queries for the randomisation module:

**A.** Random sequence of a given length

```text title="random_sequence.rinp"
--chemistry DNA
--length 12
--complement HOMO
```

`--randomise` outputs:

```text title="random_sequence.binp"
--sequence DG, DG, DG, DC, DT, DG, DA, DA, DC, DC, DC, DG
--complement HOMO
```

**B.** Prepend a given sequence with a chemistry

TODO: take this out. Prepending a sequence, nobody every used it

```text title="prepend_sequence.rinp"
--chemistry DNA
--sequence A, G, G, C, T, A, A, G, T, T, A, C
--complement HOMO
```

`--randomise` outputs:

```text title="prepend_sequence.binp"
--sequence DA, DG, DG, DC, DT, DA, DA, DG, DT, DT, DA, DC
--complement HOMO
```

## Sequence

`--sequence`: A mutually exclusive flag with [`--length`](#length). If a given sequence is given, then the length of it is already determined by the user.
The sequence is a list of comma-separated values, consisting of one of five currently implemented nucleobases: `A, C, G, T, U`.

## Length

`--length`: A mutually exclusive flag with [`--sequence`](#sequence). This requires an **integer** argument to be passed into.
It returns a randomised sequence, for a given chemistry, of the leading strand with a set length.

## Chemistry

`--chemistry`: Prompt a **chemistry type**, from a list, to query Ducque for a correctly typed sequence.
The type of chemistry is related to the **keys** of the entries given in
[Adding the backbone atoms](../user-implementation/library-tables.md#adding-the-backbone-atoms-by-which-ducque-needs-to-build).

## Complement

`--complement`: The `--complement` argument can be filled in at the randomisation stage, but its sole function is to pass the argument to the build module later.
This is an optional argument.
