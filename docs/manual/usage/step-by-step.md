# Step-by-step instructions

## How to start implementing nucleosides

This section contains a checklist of the preparation and the execution to implementing new nucleosides.

1. Perform Conformational Sampling and curate the conformation of choice to implement and build with.[^mattelaer2021]
2. Use the [xyz_pdb](xyz-to-pdb.md) function to convert the generated `xyz` to `pdb` files.
   Be mindful of the **residue** name that is given to the nucleoside, as this is important for later.
3. Appending data to the TABLES
    1. [Adding a nucleoside-linker combination](../user-implementation/library-tables.md#adding-a-nucleoside-linker-combination-to-form-a-defined-nucleotide)
       requires the user to make **residue name - inputfile.json** pairs of the nucleoside the user wants to incorporate.
    2. [Adding one or more conformations](../user-implementation/library-tables.md#adding-one-or-more-conformations-of-a-nucleoside-chemistry)
       requires the user to implement the same information as in the previous item,
       but allows the user to implement more than one different conformation of the same nucleoside.
    3. [Adding the backbone atoms](../user-implementation/library-tables.md#adding-the-backbone-atoms-by-which-ducque-needs-to-build)
       is arguably the most important aspect. This calls for the user to make **chemistry - backbone atoms** pairs.
       This query is what allows Ducque to build with any type of backbone.
    4. [Adding transmutation specifications](../user-implementation/library-tables.md#adding-transmutation-specifications)
       is only necessary to specify the full name of the chemistry for the `json` file that will be transmuted for.
    5. [Adding a nucleoside-linker pair explicitly](../user-implementation/library-tables.md#adding-a-nucleoside-linker-pair-explicitly)
       requires the user to couple the **chemistry nucleoside - chemistry linker** together.
4. Start working towards [transmuting](transmute-nucleoside.md).
   Optimise **bondangle** and **dihedral** angle values ([Defining the correct torsion angle](../user-implementation/torsion-angles.md)).
5. [Randomise](randomise.md) or query a specific sequence.
6. [Build](build.md) the desired model duplex.

## How to start implementing linkers

This section contains a checklist of the preparation and the execution to implementing new linkers.
Conformational Sampling is an option to derive the behaviour of the linker, but is by no means necessary for the model building itself.

1. Design a linker through a preferred Molecular Viewing Software (UCSF Chimera, PyMol …).
2. Designing a new linker always entails making a new nucleoside chemistry as well.

   All molecules are annotated with a specific residue name (e.g. DNA Adenosine is DA, see [`--residue`](xyz-to-pdb.md#residue)).
   A linker, by itself, needs to be associated with the nucleoside it is appended to.
   This means that the linker takes on the **residue name** of the nucleoside, to form a full nucleotide.
   Only the nucleoside `json` inputfiles have a residue name, which makes the linker much more modular to implement for other chemistries.
   The other side of the coin is that when a new linker is desired for implementation, a nucleoside chemistry has to be designed as well.

   For example: when implementing the thiophosphate linker ([Available input by default](available-inputs.md)),
   a copy of the standard DNA nucleosides were made, where the residue name was altered to be unique for the thiophosphates.
   This is because when the user queries a specific sequence, the user does this by querying specific nucleotides
   (nucleoside-linker pairs, see
   [Adding a nucleoside-linker combination](../user-implementation/library-tables.md#adding-a-nucleoside-linker-combination-to-form-a-defined-nucleotide)).
   This means that an additional nucleoside `json` needs to be created to make a nucleotide,
   since we need to be able to distinguish this type of information explicitly to the library.

   Additionally, because the thiophosphate linker has a stereocenter, two copies were made of the DNA chemistry.
   The residue name for the R-isomer became [1DA, 1DG, 1DC, 1DT], for the S-isomer it became [2DA, 2DG, 2DC, 2DT].
   This naming is arbitrary and is there to allow the distinction between regular phosphate linkers and the two thiophosphate stereomers.
   While the nucleosides [DA, 1DA, 2DA] are identical in properties, they differ by **residue name** only and the specific **linker** to which they are coupled in the library
   ([Adding a nucleoside-linker combination](../user-implementation/library-tables.md#adding-a-nucleoside-linker-combination-to-form-a-defined-nucleotide),
   [Adding a nucleoside-linker pair explicitly](../user-implementation/library-tables.md#adding-a-nucleoside-linker-pair-explicitly)).
3. Appending data to the TABLES
    1. [Adding a nucleoside-linker combination](../user-implementation/library-tables.md#adding-a-nucleoside-linker-combination-to-form-a-defined-nucleotide)
       requires the user to make **nucleoside-linker** (`json` format) pairs of the nucleoside the user wants to incorporate.
    2. [Adding one or more conformations](../user-implementation/library-tables.md#adding-one-or-more-conformations-of-a-nucleoside-chemistry)
       requires the user to implement the same information as in the previous item,
       but allows the user to implement more than one different conformation of only the same nucleoside.
    3. [Adding the backbone atoms](../user-implementation/library-tables.md#adding-the-backbone-atoms-by-which-ducque-needs-to-build)
       is arguably the most important aspect. This calls for the user to make **chemistry - backbone atoms** pairs.
       This query is what allows Ducque to build with any type of backbone. This is for both the **nucleoside** and **linker** chemistries.
    4. [Adding transmutation specifications](../user-implementation/library-tables.md#adding-transmutation-specifications)
       is only necessary to specify the full name of the chemistry for the `json` file that will be transmuted for.
    5. [Adding a nucleoside-linker pair explicitly](../user-implementation/library-tables.md#adding-a-nucleoside-linker-pair-explicitly)
       requires the user to couple the **chemistry nucleoside - chemistry linker** together.
4. Start working towards [transmuting](transmute-linker.md).
   Optimise **bondangle** and **dihedral** angle values ([Defining the correct torsion angle](../user-implementation/torsion-angles.md)).
   Generally speaking, the linker just needs to not clash with the surrounding atoms.
   The angle values here may well be an arbitrary value that evades clipping of the moieties.
   Subsequent minimisation and MD simulations will take care of the rest.
   To highlight, the phosphate and thiophosphate linkers have the exact same angle values.
5. [Randomise](randomise.md) or query a specific sequence.
6. [Build](build.md) the desired model duplex.

[^mattelaer2021]: Mattelaer CA, Mattelaer HP, Rihon J, Froeyen M, Lescrinier E.
    Efficient and Accurate Potential Energy Surfaces of Puckering in Sugar-Modified Nucleosides.
    *J. Chem. Theory Comput.* 2021, 17, 3814–3823. DOI: [10.1021/acs.jctc.1c00270](https://doi.org/10.1021/acs.jctc.1c00270)
