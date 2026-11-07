# Defining the correct torsion angle for the given nucleoside

![Workflow for defining backbone torsion angles: PES, scatterplot, uncapped nucleoside, embedded build module and iterative optimisation](figures/how-to-build-example.svg)

***A.** Potential Energy Surface (PES) of DNA Adenosine.
**B.** A scatterplot that was used to create the PES, by interpolating gridpoints.
**C.** Uncapped DNA Adenosine, ready for conversion to a `json` file.
**D.** The build module, embedded into the transmute module. This makes for easy access to building new structures after modifying angle values.
**E.** Closer look to how one can iteratively optimise backbone angle values. The complementary strand is not shown here to reduce visual noise.
**F.** Leading strand in coral, complementary strand in cornflower.
This exemplifies how the differences in backbone angle values can have a dramatic change on the resulting structure.*

In panel C of the figure above, the user can curate the Conformational Sampling[^mattelaer2021] landscape
for a specific conformation to be built into Ducque's repository.
Here is a depiction of the DNA Adenosine nucleoside ($^{2'}E$), prepared as an uncapped nucleoside, ready to be transmuted by Ducque.

In panel D, the embedded build module is shown, which is a callback to step 5 of the
[Transmute Nucleoside GUI](../usage/gui.md#transmute-nucleoside).
Here, the user loads in a finished `binp` to the **Import build file** button.
Whenever the user makes changes to the produced `json` library file, upon transmuting, the user just has to click the **Build!** button to produce a new structure,
which will build a structure with the recently modified backbone angle values.
This makes it possible for the user to iteratively adjust and build structures without having to leave the transmute GUI.

The author suggests the following protocol, as devising the correct backbone dihedrals can be quite tricky when trying out for the first time.

- It is advisable that, whenever curating conformers from the landscape (panel B), the same conformer is picked for all nucleobase types.
  It has been noticed that purines and pyrimidines tend to have similar backbone angle values respectively.
  This makes the model building simpler; the goal to a good structure is building a model that is fit for upcoming simulations.

  Additionally, it is mandated to flatten the nucleobase and glycosidic improper torsion by means of geometry optimisation at a very cheap level (i.e. HF-3c).
  Optimising at this level will cause the substituents to lie in the plane of the nucleobase.
  The torsion, however, requires to be restrained at an angle of (-179.5°).
  For example: for Adenosine, this would be the set `[C1', C8, N9, C4]`.
- As shown in panel E, the structure can start from a starting value, here arbitrarily chosen to be $\epsilon_{120^\circ}$.
  By opening your preferred Molecular Viewing Software (UCSF Chimera, PyMol …), one can select a dihedral and rotate around it,
  up to the point that the basepairs stack well with each other.
  Emulating, or at least trying to emulate, stacking conditions is a good indicator of going into the right direction.
- The leading strand is entirely dependent on the two following items:
    1. The conformation of the nucleoside prompted by the user. This is the most obvious one, but is mentioned nevertheless.
       The specific conformation the user employs in the model building will be of key importance in how the duplex turns out.
       This also means the geometry of the backbone of that conformer is important.
    2. The **leading strand** is built according to the β angle and the last two angles of the backbone, not accounting for the glycosidic bond angle.
       For example: DNA will be built along the highlighted angles: [α, **β**, γ, δ, **ε**, **ζ**, χ].

  → This means that all other angle values, except for α, are already defined with the prompted conformation.
  They are required to be filled in by the user, when transmuting, for the sake of completeness!
- The **α** and **ζ** are more important for the complementary strand's associated fitting,
  but nucleobase positioning is the prime contributor for the complementary strand build process.
  **This is why it is very important to keep nomenclature consistent in the nucleobases**[^iupac1983],
  as basepairing rests on the shoulders of correctly named atoms to infer basepairs!
- Another suggestion is to use **poly-N sequences** to build the structure.
  This will make it far easier to modify angles of one type of nucleoside at one time.

  **NB**: in order to build a **poly-N structure**, the **complement of N** needs to be implemented too!
  Though the complement does not need to have its angles optimised just yet,
  as basepairing conditions are the driving factor in fitting the complementary strand when building the molecular structure.

[^mattelaer2021]: Mattelaer CA, Mattelaer HP, Rihon J, Froeyen M, Lescrinier E.
    Efficient and Accurate Potential Energy Surfaces of Puckering in Sugar-Modified Nucleosides.
    *J. Chem. Theory Comput.* 2021, 17, 3814–3823. DOI: [10.1021/acs.jctc.1c00270](https://doi.org/10.1021/acs.jctc.1c00270)

[^iupac1983]: IUPAC-IUB Joint Commission on Biochemical Nomenclature (JCBN).
    Abbreviations and Symbols for the Description of Conformations of Polynucleotide Chains.
    *Eur. J. Biochem.* 1983, 131, 9–15. DOI: [10.1111/j.1432-1033.1983.tb07225.x](https://doi.org/10.1111/j.1432-1033.1983.tb07225.x)
