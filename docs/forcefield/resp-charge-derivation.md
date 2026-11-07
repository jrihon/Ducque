# RESP Charge derivation

## The R.E.D. server

The R.E.D. server is a tool to automate the RESP charge derivation.
However you'd want to derive RESP charges is up to you.
One can easily do it on their local system and doing it manually gives you a very good understanding of it all
(see [RESP charge derivation on your local system](#resp-charge-derivation-on-your-local-system)).

The R.E.D. server will run geometry optimisations of your molecule before calculating the atomic charges.
This optimisation is by far the most expensive step of it all.
It is advisable to add an extensive `System.config` and a `Project.config` to your inputs, to optimise the accuracy of your calculations.

The R.E.D. server does not need registration.

Computational chemistry available through the following packages:

- Gaussian (Advisable to stick to this software! Commonly used in the AMBER community)
- GAMESS-US
- Firefly
- RESP
- [R.E.D Server Development news](https://upjv.q4md-forcefieldtools.org/REDServer-Development/news.php)

### R.E.D. server documentation

Read these publications as they contain quite some information.

- Documentation of the RESP charge derivation:
  *A Well-Behaved Electrostatic Potential Based Method Using Charge Restraints for Deriving Atomic Charges: The RESP Model* (Bayly et al.)[^bayly1993]
- Documentation of the R.E.D. server:
  *R.E.D. Server: a web service for deriving RESP and ESP charges and building force field libraries for new molecules and molecular fragments*
  (Vanquelef et al.)[^vanquelef2011]
- Improvements on the R.E.D. server:
  *The R.E.D. tools: advances in RESP and ESP charge derivation and force field library building* (Dupradeau et al.)[^dupradeau2010]
- PyRED stable interface: [PyRED server](https://upjv.q4md-forcefieldtools.org/REDServer-Development/)
- R.E.D. Server homepage: [Homepage Q4MD-Forcefieldtools](https://upjv.q4md-forcefieldtools.org/)
- Here you can find how to prepare input files for a first R.E.D. Server development:
  [Inputfile preparation HowTo](https://upjv.q4md-forcefieldtools.org/Tutorial/Mini-HowTo-InputFiles.pdf)

### Basics

The tutorials of the R.E.D. server contain a section on how the RESP derivation scheme works.
This can help one with understanding how restraints and equivalencing works, either through the server or on the local system.

Be sure to check out the **Tutorial: [Updated Tutorial](https://upjv.q4md-forcefieldtools.org/Tutorial/Tutorial-4.php)!**

### Configuration files

The following files are required to run a R.E.D. server job:

- System.config: [Link to System.config](https://upjv.q4md-forcefieldtools.org/REDServer-Development/Documentation/System.config)
- Project.config: [Link to Project.config](https://upjv.q4md-forcefieldtools.org/REDServer-Development/Documentation/Project.config)
- Some added information on the format of the input files: [README HowTo](https://upjv.q4md-forcefieldtools.org/REDServer-Development/Documentation/readme.txt)

Everything here should guide you to be able to manage and generate a very simple charge derivation that is compatible with AMBER, GROMACS, OPLS or CHARMM.

## RESP charge derivation on your local system

### Small set-up and reasoning behind the RESP calculation

The ORCA manual is available on their repository over at the [ORCA Forum](https://orcaforum.kofo.mpg.de/app.php/portal).
You will have to make a free account and then you get access to the Downloads page and there are several folders.
Download the latest version of ORCA, together with the respective version of Open MPI, also free of charge ([Open MPI official homepage](https://www.open-mpi.org/)).
Installing ORCA is easy as you will only download the binaries, no need to compile anything by hand (looking at you, GAMESS).

Go to the *ORCA Manual* folder to download the correct manual; you will need it. It contains the installation procedure!

Furthermore, the charge derivation will be done by ORCA in a first step, to compute for *HF/6-31G\**-derived wave functions.
This will garner us the data to plot the Electrostatic Potential (ESP) on a grid that is generated through the MSMS software,
see [Generation of the grid and plotting ESP on the grid](#generation-of-the-grid-and-plotting-esp-on-the-grid).
The ESP-loaded grid will be fed to the RESP (`$AMBERHOME/bin/resp`) script, provided with AMBER.

The **`ff/scripts/`** directory contains multiple scripts and files:

- `msms_setup.sh`: To help and direct the extraction of the MSMS software, needed for the MK scheme.
- `reorient_molecule.py`: Used to do a reorientation of the molecules to produce more robust charges for resp.
- `template_Orca.inp`: a template ORCA input file for the *HF/6-31G\** computation.
- `OrcaExec.sh`: A shell script to run ORCA commands. It is more practical than writing out full commands.
  On top of that, the shell command `orca` exists on Unix systems and it avoids confusion by calling the correct script.
- `MK_PopAnalScheme.py`: Used to generate the ESP-loaded grid for RESP.

To first calculate the charges, we need the correct level of theory and basis set.
For consistency reasons, we will have to assess the single point energies through a Hartree-Fock level, with an older basis set.
This is done because the AMBER forcefield has, in the nineties, been parametrised using this. More specifically with the *6-31G\** basis set.
At that time, they saw the most consistent results happen with this computation and stuck to it, so newer AMBER forcefields have to be designed around this.
Luckily for us, calculations at the HF-level only take up several minutes nowadays.

Fun fact: the *HF/6-31G\** computation happens in gas phase, while we usually do our simulations in a solution.
Conveniently, the error on the single point energy with this computation is about the same as the difference in charge distribution between gas and solution phase,
so it cancels out the error!
[Better explanation on the use of this level of theory and basis set.](http://ambermd.org/tutorials/advanced/tutorial1/section1.htm)

How RESP works is that the ESP-loaded grids are supplied to the compiled `resp.F` script,
along with added information on possible restraints and equivalencing on the molecules.
This is especially important when we work with fragmented molecules, like amino acids and nucleic acids.
See [Restrained Electrostatic Potential (RESP) with AMBER](#restrained-electrostatic-potential-resp-with-amber) for a small introduction to this.
For a complete analysis of how a nucleoside forcefield is created, we refer to the following publications:
*RESP methodology: Charge derivation for DNA, RNA and proteins*[^cieplak1995] and *AMBER FF for modified RNA*[^aduri2007].

In 1995, a second generation of AMBER forcefield parameters came to be in existence (highly advisable read!): *2nd generation AMBER FF*[^cornell1995].
There is also the [RESP Frequently Asked Questions (FAQ) page](http://ambermd.org/Questions/resp.html).

### Electrostatic Potential (ESP) calculation with ORCA

We need to have the information on orbital densities and molecule coordinates in order to calculate for the ESP charges,
which have been mapped on a specific grid (more on this later).
This in turn will be used to derive RESP charges through the `resp` program (part of AmberTools).

The **`ff/scripts/`** directory contains a script named **`MK_PopAnalScheme.py`** that does all the following for you.

Whenever we do a Single Point energy evaluation we generate a set of files with information on densities, orbitals,
coordinates (xyz) of the molecules' atoms and such. These files are required by the `orca_vpot` program.

What follows is the output of `orca_vpot` whenever no arguments have been prompted.

```text
   Interactive use:
   $ orca_vpot  GBWName PName XYZName POTName [BaseName of DensityContainer]

       GBWName  = GBW File that contains the orbitals/coordinates/basis
       PName    = Name of density (or density container) (must match the GBW file basis set!)
       XYZName  = File that contains the grid coordinates to evaluate V(r) for
       POTName  = Output file with V(r)
       DensName = BaseName of DensityContainer (if different from BaseName of GBWFile)
```

### Generation of the grid and plotting ESP on the grid

Since ORCA has not (yet) implemented the Merz-Kollman population analysis scheme, we are doomed to do it ourselves.
Fear not, for yours truly has scripted and taken the brainwork out of it. Download the tar file, pop a button and you're done.
A shell script will be provided to get the **MSMS** program going. You will have to specify the name of the executable (binary) for the script to work.

- [Molecular Graphics Laboratory (MGL) Homepage; The Scripps Research Institute](https://ccsb.scripps.edu/mgltools/#msms)
- [MSMS software download page](https://ccsb.scripps.edu/msms/downloads/) → Download the linux tar file!
- [MSMS manual page UPDATED](https://ccsb.scripps.edu/msms/documentation/) → The manual page gives extra context and guidance on how to use the software.
- [MSMS manual page DEPRECATED](https://www.scripps.edu/sanner/html/msms_man.html)

There are multiple ways to do this, but I did it like this.
Retrieve the `msms_setup.sh` script from the `ff/scripts/` directory or copy the following code into a script (i.e. `msms_setup.sh`)
and run the script with the following line:

```bash
bash msms_setup.sh
```

```bash
#!/bin/bash

echo "Change the name of the LinuxX_X.X.X to the correct tar file name."
echo "Then remove or comment out the exit command from the script when succesful."
exit
mkdir $HOME/msms

tar xvzf $HOME/Downloads/msms_i86_64LinuxX_X.X.X.tar.gz --directory $HOME/msms
cd $HOME/msms

cp msms.x86_64LinuxX.X.X.X msms
```

Add the following line to your `$HOME/.bashrc` file manually.

```bash
export PATH=$PATH:$HOME/msms/
```

The workflow is as follows:

1. Reorient the molecules (onto a predefined plane) using `reorient_molecule.py`.
   Requires argument *nucleoside* or *linker*. Read its docstring for more information.
2. Carry out a *HF/6-31G\** computation on the molecule.
   There is a `template_Orca.inp` file present in **`ff/scripts/`**, this should be altered to your specific inputs.
   This directory also contains a `OrcaExec.sh` script to run ORCA with the correct flags.
   It saves a bit of syntax hassle and disarms ambiguous calls (looking at you `/usr/bin/orca`).
3. Run the `MK_PopAnalScheme.py` script with the correct basename of the ORCA `*.out` file as an argument.
   The docstring, in the script, details how it works and what the commandline-arguments should be.

**NOTE:** Make sure everything is installed as detailed! The script is mainly hardcoded, with respect to runtimepaths, so small errors are quickly made.
It should generally not be difficult to figure out arising errors as they are more often than not runtimepath-errors.
Everything has been well documented for your convenience.

### Restrained Electrostatic Potential (RESP) with AMBER

Now that we have the ESP charges mapped onto the grid
(by using `orca_vpot` to plot the *HF/6-31G\** calculated wave functions on our own generated grid), we can use this as inputs for the `resp` program.

The `resp` program essentially uses a least-squared-fitting algorithm to derive the RESP charges
(which are just a specific type of atomic charges) from the ESP-loaded grid generated earlier.
See the Ducque publication on a graphical representation of this process.

- [FAQ about resp](http://ambermd.org/Questions/resp.html)
- [Advanced tutorial on resp](http://ambermd.org/tutorials/advanced/tutorial1/section1.htm)

The derivation of RESP charges takes some getting used to but essentially boils down to this:

1. Free derivation of atomic charges, unless otherwise restrained with predefined charges from a library.
2. Restraining and equivalencing of degenerate hydrogens (in methyls and methylenes), by using the charges from the first RESP derivation.

Here are a number of **CAVEATS** I had when generating charges for the nucleosides:

1. Make sure to format exactly as you see in the tutorial. The `resp` script was written in Fortran and parses on a column-basis.
   I have included a big docstring-like document at the end of this manual (see [RESP DocString](miscellaneous.md#resp-docstring)).
   Especially the information about formatting is in there. I translated the Fortran formatting to python formatting.
    - Most columns (section 5, 6, 7.1, 8.1) are formatted as `"%5d"`.
      Which means every integer you input only gets a maximum of 5 characters (length) total.
      This does not mean five whitespaces and then an integer. This means your value should fit in the span of five characters.
    - There should be left an empty line as the last line of the file.
    - There should be an empty line after every atom's information (section 3, 4, 5, 6).
    - Between the molecules' information and the restraining of atoms to a certain charge should be a blank line (Between Section 6 and 7 or 8).
      Section 7 and 8 constitute the same type of information and will be considered as similar sections.
    - Between restraining of charges (section 7 and/or 8) and inter-molecule equivalencing should be a blank line.
    - Indent the first column, of section 7.1/8.1, with four spaces.
    - If no charge restraining is required (so no section 7 and 8), then leave two blank lines between molecule information (section 6)
      and the equivalencing part (section 9).
2. The output of `resp` is badly written. I have spent hours trying to analyse its output and make sense of it.
   For example it takes the name of the last molecule prompted and will use that title as the only title included in the output of the file at the punch (`-p`) output.
   The output (`-o`) on the information of the charge restraining is also not clear.
   That information is outputted together with information of the last molecule prompted (section 6).
   This was especially confusing, as it gives the impression the charge restraining was not prompted correctly.
3. A good check-up of having a successful derivation of the atomic charges is summing the charges, produced by `resp`,
   of the restraint process (section 7 or 8) and see if they add up to the ultimate charge they were fitted to.
   For example, neutral DNA fragments add up to a net charge of zero (0).
   Terminal 5' and 3' fragments are added up together to a charge of minus one (-1). Internal fragments add up to a net charge of minus one (-1).
4. Definitely recommended, to see if the charge derivation and the MK scheme were successful,
   is to compare the results from your output to that of the following file: `$AMBERHOME/dat/leap/prep/nucleic10.in`
   and see if the charges are in the same order of magnitude.
   Important is that they should not specifically be positive or negative. Some carbons can be very slightly positive, some a bit negative.
   It all depends on the nucleic acid chemistry you are characterising. Just try to make sense of the charges in general.

The conversion from `resp` to the `.prep` file format is automatable.
Since it is prone to errors and, unless you are asked a bunch of set input parameters,
it will be faster done manually than trying to figure out which input parameters are correct.
There is a small script in `ff/scripts/` that parses all the charges per prompted molecule to a separate file,
so you can more easily mix and match the required charges of the respective atoms to make your nucleo(t)(s)ide fragments.
Advisable to use the `cat` command to concatenate all the different files into one big `.prep` file.

Amino acids and nucleic acids are divisible into fragments, to represent them well depending on their manner of appearing.
For all nucleotides specifically, there will be a set of four charges per respective nucleobase.

- A 5' (OH) with a 3' (O⁻) fragment; 5-prime ending fragment (DX5)
- A 5' (O⁻) with a 3' (OH) fragment; 3-prime ending fragment (DX3)
- A 5' (O⁻) with a 3' (O⁻) fragment; A fragment (a.k.a residue) in the nucleic acid strand (DX)
- A 5' (OH) with a 3' (OH) fragment; A neutral fragment, not incorporated in a nucleic acid strand (DXN)

Very important, when parsing the charges, is that when combining charges of the linker and the nucleoside,
we only include the charges of the atoms that have not been restrained to a set charge.
This means that for the 3'-oxygen in a standard nucleotide, one uses the charge of the 3'-oxygen that has been computed for in the dimethylphosphate linker.

In the case of the morpholino sugar moiety, for the nitrogen in the ring,
we use the charge computed for in the N,N-dimethylaminophosphoramidate linker moiety,
since this nitrogen was not restrained with a group of atoms to the net charge of zero (0).
More on this in the articles on the R.E.D. server ([R.E.D. server documentation](#red-server-documentation))
and the 2nd generation forcefield ([Small set-up and reasoning behind the RESP calculation](#small-set-up-and-reasoning-behind-the-resp-calculation)).

[^bayly1993]: Bayly CI, Cieplak P, Cornell W, Kollman PA.
    A well-behaved electrostatic potential based method using charge restraints for deriving atomic charges: the RESP model.
    *J. Phys. Chem.* 1993, 97, 10269–10280. DOI: [10.1021/j100142a004](https://doi.org/10.1021/j100142a004)

[^vanquelef2011]: Vanquelef E, Simon S, Marquant G, Garcia E, Klimerak G, Delepine JC, Cieplak P, Dupradeau FY.
    R.E.D. Server: a web service for deriving RESP and ESP charges and building force field libraries for new molecules and molecular fragments.
    *Nucleic Acids Res.* 2011, 39, W511–W517. DOI: [10.1093/nar/gkr288](https://doi.org/10.1093/nar/gkr288)

[^dupradeau2010]: Dupradeau FY, Pigache A, Zaffran T, Savineau C, Lelong R, Grivel N, Lelong D, Rosanski W, Cieplak P.
    The R.E.D. tools: advances in RESP and ESP charge derivation and force field library building.
    *Phys. Chem. Chem. Phys.* 2010, 12, 7821–7839. DOI: [10.1039/c0cp00111b](https://doi.org/10.1039/c0cp00111b)

[^cieplak1995]: Cieplak P, Cornell WD, Bayly C, Kollman PA.
    Application of the multimolecule and multiconformational RESP methodology to biopolymers: Charge derivation for DNA, RNA, and proteins.
    *J. Comput. Chem.* 1995, 16, 1357–1377. DOI: [10.1002/jcc.540161106](https://doi.org/10.1002/jcc.540161106)

[^aduri2007]: Aduri R, Psciuk BT, Saro P, Taniga H, Schlegel HB, SantaLucia J.
    AMBER Force Field Parameters for the Naturally Occurring Modified Nucleosides in RNA.
    *J. Chem. Theory Comput.* 2007, 3, 1464–1475. DOI: [10.1021/ct600329w](https://doi.org/10.1021/ct600329w)

[^cornell1995]: Cornell WD, Cieplak P, Bayly CI, Gould IR, Merz KM, Ferguson DM, Spellmeyer DC, Fox T, Caldwell JW, Kollman PA.
    A Second Generation Force Field for the Simulation of Proteins, Nucleic Acids, and Organic Molecules.
    *J. Am. Chem. Soc.* 1995, 117, 5179–5197. DOI: [10.1021/ja00124a002](https://doi.org/10.1021/ja00124a002)
