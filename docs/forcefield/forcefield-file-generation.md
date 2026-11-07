# Forcefield file generation

The section in the AMBER manual titled ***Reading and modifying Amber parameter files*** (under System preparation)
gives a detailed explanation on parameter topology (`.prmtop`) and forcefield modification (`frcmod.` or `.dat`).
For the sake of being complete in the manual, I will give a brief overview of their importance.
When one wants to do an extensive parametrisation of multiple conformations, it is advised to also follow this tutorial (after generating RESP charges):
[ParamFit Tutorial](https://ambermd.org/tutorials/advanced/tutorial23/index.php).
This program will generate a more suitable `frcmod.` file to better estimate the conformations' potential in molecular mechanics.

The `$AMBERHOME/dat/leap` directory contains three directories of interest.

- `cmd`: Contains the leaprc (leap-run-command) files that are sourced to import forcefields into tLEaP when generating a set of parameters.
  These leaprc files in their turn call upon `.prep`, `frcmod.`, `.lib`, `.dat` files.
- `parm`: Contains the `frcmod.` files that contain information of specific torsion angles.
  This contains information of bond stretching, angles bending and phase torsion.
  The `.dat` and `.lib` file formats are a variation on the `frcmod.` file, but detail the same parameters.
- `prep`: Contains the `.prep` files that contain information about connectivity,
  atom typing (required for atom geometry and orbital hybridisation) and charges.
- Here is a link to the standard file formats of AMBER. Note the "Main parameter set" and "Parameter modification file" are of special interest here.
  [AMBER File Formats](https://ambermd.org/FileFormats.php)
- The following publication has added information on how the original AMBER FF were built,
  which includes fundamental information on the given parameters in the `.prep` and `frcmod.` files:
  *A new force field for molecular mechanical simulation of nucleic acids and proteins*[^weiner1984].
  It is an older publication for sure, but it checks out.

Let us say we start off tLEaP.
The first thing that is instantiated is the path to the different directories are set to the runtimepath of your environment,
meaning you can source (or call) any files within those directories without having to explicitly type the path to every forcefield file you need.
For the sake of simplicity, I assume the `$AMBERHOME` has been installed inside the `$HOME` directory
(the real output will never show exported PATHs, but the absolute path to the files and directories. I am generalising here):

```text
-I: Adding $HOME/amber18/dat/leap/prep to search path.
-I: Adding $HOME/amber18/dat/leap/lib to search path.
-I: Adding $HOME/amber18/dat/leap/parm to search path.
-I: Adding $HOME/amber18/dat/leap/cmd to search path.

Welcome to LEaP!
(no leaprc in search path)
> source leaprc.DNA.OL15
----- Source: $HOME/amber18/dat/leap/cmd/leaprc.DNA.OL15
----- Source of $HOME/amber18/dat/leap/cmd/leaprc.DNA.OL15 done
Log file: ./leap.log
Loading library: $HOME/amber18/dat/leap/lib/DNA.OL15.lib
Loading parameters: $HOME/amber18/dat/leap/parm/parm10.dat
Reading title:
PARM99 + frcmod.ff99SB + frcmod.parmbsc0 + OL3 for RNA
Loading parameters: $HOME/amber18/dat/leap/parm/frcmod.DNA.OL15
Reading force field modification type file (frcmod)
Reading title:
OL15 force field for DNA (99bsc0-betaOL1-eps-zetaOL1-chiOL4) see http://ffol.upol.cz
>
```

What happens here is that we just read into tLEaP's environment the forcefield files that contain all the information on DNA parameters.
This information entails bonded and non-bonded term interactions.

The leaprc (analogous to the Linux system `~/.bashrc`) calls upon the respective `.prep` and `frcmod.` files
and also lets tLEaP know that sometimes, conventional atom names are switched up here and there.
In order to recognise them, it will assign synonyms. Like the oxygens in the phosphate backbone can be named either "OP1, OP2" or "O1P, O2P".
To keep everything in check, and because conventions change every so often (thereby defeating the purpose of a convention haha),
we want consistency in our parameters.

Next to this, the `.lib`, `frcmod.` and the `.prep` are read by tLEaP to gather all the necessary information.
The reason I did not include the `.lib` file format in the explanation is because I do not think it is a relevant format to create current-day forcefields with,
but mainly with the latter two formats. You can be just fine creating and employing a forcefield without the `.lib` format files.

## Parameter Topology .prmtop

Required to run your simulations with. These are always accompanied with a `.crd` or `.mdcrd` (coordinate) file.
The prmtop points to coordinates, that represent the atoms in the system, and assigns charges, bondlength, angle, 1-4 torsion parameters, atom types and so forth.
A more detailed depiction of what entails a `.prmtop` and a `.crd` file is here:
[AmberMD tutorial: Fundamentals of LEaP](https://ambermd.org/tutorials/pengfei/index.php).

The `.prmtop` and the `.crd` files are outputs of the tLEaP program,
which are convenient formats for the *sander* and *pmemd* simulation engines to read all-atomic data from.
These file formats are created from a set of different files and are most commonly derived from `.frcmod` and `.prep` files,
which are inputted into tLEaP to generate parameters for specific models.

## Preparatory files .prep

The `.prep` file contains connectivity information (which are the internal coordinates) and the charges of the molecules.
This file also contains the information on residue names as well, which is one of the ways in which tLEaP recognises molecules to parametrise.
Lastly, this file format contains the AMBER filetree structure, which is AMBER's way of making a graph-like structure out of a molecule.
See [Tree Structure: Topological Types](prep-file.md#tree-structure-topological-types).

Additionally, it is also possible to add loop closing statements (for cyclic-like moieties) and to impose a set of dihedral restraints through the improper flag.

One can generate a standard `.prep` file through the `antechamber` script provided with AMBER,
though for complex molecules like nucleic acids the antechamber script will not suffice.

Instead this part will have to be done manually.
I have supplied the repository (`ff/scripts/`) with a small script that helps with formatting of the prep file, named `format_prep.py`.
It seems that antechamber natively assumes the tree structure of some molecules and it messes up the `.prep` file this way.
It is also not a hundred percent correct on assigning atom types. Therefore, it will have to be done manually.

The good thing about having a model builder, like Ducque, at hand, is that there is a lot of information in the `.prep` file that can be skipped.
For instance, the connectivities, the values for the bonds, angles and dihedrals can be skipped.
These get overridden by a prompted model and the AMBER Tree Structure.

The important columns are the atom names, amber tree structure, amber/gaff atom typing and of course the charges column.

### Practicals

I have a directory with a total of sixteen (16) directories, containing a neutral fragment, a head fragment,
a tail fragment and an internal fragment for all (4) nucleosides.
See [Restrained Electrostatic Potential (RESP) with AMBER](resp-charge-derivation.md#restrained-electrostatic-potential-resp-with-amber), right at the end.

The script `format_prep.py` requires two arguments. The first is a file consisting of four columns, like here below.
The second argument is the residue name of the molecule.
For the people comfortable in python, you can mix and match however you read in the file or how the input file is formatted,
I just wrote something that works for me as a comfortable structure of inputs.

```text
ATOM NAME - TREE STRUCTURE - ATOM TYPES - CHARGES
HO6'    M    HO    0.446871
 O6'    M    OH   -0.633980
 C6'    M    CI   -0.032336
H6'1    E    H1    0.087206
H6'2    E    H1    0.087206
 C5'    M    CT    0.128480
 H5'    E    H1    0.110299
 O5'    S    OS   -0.273606
...    ..    ..      ....
 C2'    B    CT   -0.162348
H2'1    E    H1    0.158226
H2'2    E    H1    0.158226
HN3'    M    HN    0.364026
```

## Forcefield modification frcmod

`$AMBERHOME/dat/leap/parm` is a directory that contains all utilised standard FF parameter files for the AMBER forcefields.
The `frcmod.` files contain all information on the physical parameters of the molecule that one wants to simulate.
More specifically it assigns bond lengths, bond angles and torsion (dihedral) parameters.
Non-bonded (1-4) parameters are also defined here and these are less likely to be defined by the user.
Mainly torsion parameters are particularly difficult to get right and are most subject to tweaking according to experimental (NMR) and vibrational data.
When setting up a forcefield for a modified nucleic acid,
it is customary to get most parameters readily from these `frcmod.` files and tweak the ones not conforming with the standard DNA or RNA nucleic acid.

The `frcmod.` files contain information on energy barriers for all physical parameters.
These energy barriers are energetic values (assessed in molecular mechanics, expressed in kcal/mol)
that need to be crossed in order for a different equilibrium value to be eligible.
For example: the different puckering modes of DNA (North and South conformation) are defined as different energy barriers
and the puckering can only flip from one to another whenever the energy of the system permits the molecule to traverse the transitional path.
That is to say, some puckering modes are more likely to appear when hybridising with a certain chemistry
and that specific chemistry is likely to be more stringent or less likely to shift conformations,
thereby forcing DNA to shift and in the process the energy barriers will be crossed due to efficient interactions being made
(hybridisation and stacking interactions are strong contributors to traversing transitional states).

As this is a rather abstract concept, I suggest you read the literature on [paramfit](https://ambermd.org/tutorials/advanced/tutorial23/index.php)
to familiarise yourself with fitting of the MM parameters to QM values.
We will use Paramfit later on to evaluate certain types of bonds (like the phosphoramidate bond) ourselves,
together with the morpholino-'sugar' of our new modified nucleoside.

One can generate a standard `frcmod.` file through the `parmchk2` script provided with AMBER, typically only for small molecules.
Yet again, nucleic acids are far too complex to be handled by scripts meant for standard small molecules, so we will have to craft it ourselves.
Luckily, nucleic acids have been well-described and we can take a lot of the parameters directly from the RNA and DNA `frcmod.` files.
This is especially true for chemistries that have been 2' substituted, like 2'-O-Methyl RNA.

### Bonded Interaction Terms

The AMBER manual has an extensive documentation on the ins and outs of the `frcmod.` section (AMBER18, section 14.1.6).
This is a brief reiteration for the sake of completion.

The way the `frcmod.` file is formatted is by specific columns. A couple abbreviations are listed here to keep the list short:
Atomtype (AT), force constant (FC, this is the energy barrier earlier discussed),
equilibrium value (EV, as length / degrees / phase), divider (DIV), periodicity (PERIOD).

- **Bond Stretching**
  `AT-AT FC EV      DOCUMENTATION AND COMMENTS`
- **Angle Bending**
  `AT-AT-AT FC EV      DOCUMENTATION AND COMMENTS`
- **Torsion**
  `AT-AT-AT-AT DIV FC EV PERIOD      DOCUMENTATION AND COMMENTS`

What this essentially entails is that the simulation engines will read the `.prmtop` file,
which is a concatenation of `.prep` and `frcmod.`, and read all the charges and geometric properties in.
Whenever the simulation goes on, the structures and puckering will be evaluated and railed in whenever they go out of bounds.
Or a different energetically favourable structure can be assumed.

### Periodicity

Since this was a rather confusing definition for myself.
(source: `frcmod.*` documentation, slightly modified to update words and exclude fortran-jargon)

```text
The periodicity of the torsional barrier.
NOTE: If periodicity is less than 0.0, then the torsional potential is assumed to have
      more than one term, and the values of the rest of the terms are read from
      the next lines until a positive periodicity is encountered.
      The negative value of periodicity is used only for identifying the existence of
      the next term and only the absolute value of periodicity is kept.
```

This means that if the torsion term has more than one possibility, then the periodicity has a negative integer value.
For however many torsion angles exist for that dihedral, we will encounter a negative integer (by convention, in descending order) until a positive value is encountered.
All the torsion angles and their respective FC are read into memory.
Although a negative periodicity entails that the next line will be read too for the same dihedral set,
the absolute value of the periodicity will also be used in the equation.

In more mathematical terms, the periodicity is used to describe the torsional potential.
In other words, the equation below evaluates the [energy of the dihedral's] torsion through a Fourier series, which is expressed as a cosine function.
In some textbooks, the periodicity is sometimes referred to as the multiplicity (this confused me a whole lot, since these two are interchangeable).

As I have read, the periodicity stands for the amount of energetic minima that the torsion angle can comprise of.
In a range from [0°, 360°], following the potential of the molecule when varying the dihedral at this range,
we will see $n$ amount of minima in the potential wave function.
[sciencedirect.com - periodicity explained](https://www.sciencedirect.com/topics/chemistry/torsional-potential)
or [Ideas of Quantum Chemistry, p.287-289](https://www.sciencedirect.com/book/9780444522276/ideas-of-quantum-chemistry).

Here below, we have the equation by which the AMBER engine evaluates the potential of molecules through Molecular Mechanics.

$$
\begin{aligned}
V_{AMBER} & = \sum^{n_{bonds}}_{i} b_i(r_i - r_{i,eq})^2 \\
          & + \sum^{n_{angles}}_{i} a_i(\theta_i - \theta_{i,eq})^2 \\
          & + \sum^{n_{dihedrals}}_{i} \sum^{n_{i,max}}_{n} \frac{V_n}{2} \left[1 + \cos(n\phi_{i} - \gamma_{i,n})\right] \\
          & + \sum^{n_{atoms}}_{i<j} \left(\frac{A_{ij}}{r^{12}_{ij}} - \frac{B_{ij}}{r^{6}_{ij}}\right) \\
          & + \sum^{n_{atoms}}_{i<j} \frac{q_i q_j}{4\pi\epsilon_0 r_{ij}}
\end{aligned}
$$

The part about the dihedral torsion (third term) can be aliased as the following:

$$
\sum_{i}^{n_{atoms}} FC \left(1 + \cos(PERIOD \cdot \phi_i - EV_{eq})\right)
$$

For AMBER to evaluate the conformers properly during the simulation,
the K-term values in the BONDED ENERGY terms of the Molecular Mechanics Energy Evaluation should properly represent the energetic landscape of the molecule.
With respect to nucleic acids, this means especially proper evaluation of the puckering conformers during the simulation (see eq. 14.1 in the AMBER18 manual).
It should be noted that the $k_b$ (bond stretching), $k_{\theta}$ (angle bending) and $V_n$ (phase torsion) parameters are semantically the same concept,
but are expressed in different units. I have aliased these to the FC-value mentioned in the list earlier to avoid confusion.

## Paramfit

Paramfit has definitely gained my favour above the `mdgx` software provided in the AmberTools package, but this is more a personal preference.
The AMBER tutorials recommend `mdgx` and noted paramfit as deprecated, but we do not care for those recommendations.

The way paramfit works is that the user needs to supply QM derived single point energies,
a preliminary `.prmtop` and a set of coordinates of the conformers that have had their energy evaluated.

$$
f(N, E_{QM}, K) = \sum^{N}_{i=1}\left[(E_{MM}(i) - E_{QM}(i) + K)^2\right]
$$

More likely than not, the bond stretching and angle bending terms in the `frcmod.` need little to no re-evaluation.
Unless you are defining a new atom type, you can just reuse older bond and angle terms,
these are very well established and vary very little to not at all in unique models.
What the user needs to do is supply a `frcmod.` and/or `.dat` and/or `.lib` file that contains the bonded term interactions for your molecule of interest.
To generate a `.prmtop` file, the user only needs to supply bond stretching and angle bending terms, if some are missing.
The torsion terms will be supplied in the second step of the Paramfit process and will be evaluated afterwards.

In the torsion term, `Paramfit` is able to define the energy barrier ($K_p$), a phase ($\phi$) as well as the periodicity ($N_p$).
Of course not all torsion terms need re-evaluation! Only terms that are new and have not been readily defined yet, or need some serious tweaking.
In the example of a new type of nucleoside, where the sugar has been modified, an extensive amount of terms need to be added.
When looking at a modified nucleobase, usually the amount of terms are toned down and not at all complicated to define,
especially when just substituents have been added or modified on the nucleobase itself.

Be VERY VIGILANT HERE! Oftentimes, you might want to add several torsion terms yourself if they are missing from the preliminary `.prmtop` file.
Look up in literature if they are available. Perhaps the GAFF data files will contain some suitable enough for initial fitting.

### Gathering parameters

There are three major steps in acquiring suitable FF parameters, usually torsion terms, for your forcefield.

1. Fitting the K-value of the system.
   This K-value will fit the MM energy values, determined by the AMBER energy potential equation, to the given QM energy potential.
   This is an important, but rapid step. This will give an initial search value for the following fitting of the torsion terms.
2. The input of parameters is defined by Paramfit themselves.
   With the given `.prmtop`, the user will be prompted to fill in their wishes and will receive an exhaustive list of dihedrals (chosen by paramfit)
   that the user may or may not want to fit. This will output a list and is used for the next step.
   Note that if you have this file, you do not need to repeat it again and again, unless you require more or less torsions to be fitted.
   I would advise here to mainly fit the energy barrier ($K_p$) and the phase ($\phi$).
3. The final step will entail the actual fitting of the MM to the QM energies. This step can be repeated as often as one wants.
   You would analyse the output and based on it, you would remove a conformer from the list,
   weight them to zero (practically commenting it out) to get a fit as best as one can.
   When you plot the MM and the QM values, do not be alarmed to see the curve itself spiking more severely than before;
   the commented conformers have been exactly accounted for according to the given weight.
   Also, when one does not easily reach a great fit of the curve,
   you might want to restart all over and look for any missing torsion parameters in your initial `frcmod.` file
   or perhaps you might have missed adding one or more in the step where one defines the torsion parameters to be fitted.

The goal of Paramfit is to tweak and complete the forcefield at hand.
Do not be afraid to tweak some torsion parameters that already exist (if necessary), though be wary of which ones you actually decide to change.
You will notice an especially large size of torsion parameters needed when fitting, let's say, the parameters for a new sugar chemistry in a modified nucleic acid.
Though more often than not, one will build upon the existing DNA and RNA forcefields,
it is best to use as much as one can from those, as they have been extensively evaluated.

**CAVEAT**: If the outputted `frcmod.*` file has torsion terms where the energy barrier ($K_P$) is absurdly high
(usually one does not encounter values above 2 or 3), make sure you have added all required torsion terms to your preliminary `frcmod.` file.
If not enough bonded terms have been supplied, the AMBER engine will not be able to accurately depict the potential of the molecule,
as it is missing some crucial information.
If the problem still persists, see how good you can get the fitting done without the term containing the high $K_P$.

Here is the link to the information and tutorials on paramfit itself:

- [Paramfit Tutorial (AMBER)](https://ambermd.org/tutorials/advanced/tutorial23/index.php)
- [Paramfit Tutorial (Walker MD Lab)](http://ambermd.org/tutorials/advanced/tutorial23/index.html#Introduction)
- Paramfit publication, JCC[^betz2015]

[^weiner1984]: Weiner SJ, Kollman PA, Case DA, Singh UC, Ghio C, Alagona G, Profeta S, Weiner P.
    A new force field for molecular mechanical simulation of nucleic acids and proteins.
    *J. Am. Chem. Soc.* 1984, 106, 765–784. DOI: [10.1021/ja00315a051](https://doi.org/10.1021/ja00315a051)

[^betz2015]: Betz RM, Walker RC.
    Paramfit: Automated optimization of force field parameters for molecular dynamics simulations.
    *J. Comput. Chem.* 2015, 36, 79–87. DOI: [10.1002/jcc.23775](https://doi.org/10.1002/jcc.23775)
