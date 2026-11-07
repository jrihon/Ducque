# Transmute nucleobase (tbase)

## PDB

`--pdb`: Relative pathing to the `pdb` file of the molecule that the user wants to incorporate into the library.

Requires the full nucleobase, preferably with hydrogens. The base of the molecule needs to be stripped from the hydrogen.

## Chemistry

`--chemistry`: The **chemistry** denotes the type of modification the molecule bears.

I guess this is something like adenine, pseudo-uracil, just a string that cannot get validated.

The chemistry has direct relations to its **key** in

RE: this line above really fucking sucks haha who thought of htis.

## Moiety

`--moiety`: The Ducque software distinguishes three parts of a molecule, being the **nucleoside** and the **linker**, and now the **nucleobase**.

I guess it's just one of these three

## atoms

`--atoms`: lets user decide prompt which atoms are at the base of the molecule and form the plan to flip over

- e.g. `--atoms N9, N1, C2`
