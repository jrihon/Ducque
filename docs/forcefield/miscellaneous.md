# Miscellaneous

## AMBER atom types with GAFF alternative atom types

### Introduction

See the following link for the source:
[The AMBER Forcefield Atom Types](http://www.uoxray.uoregon.edu/local/manuals/biosym/discovery/General/Forcefields/AMBER.html#C_all)

The standard AMBER forcefield, which is attributable to Kollman and coworkers (Weiner et al. 1984[^weiner1984], 1986[^weiner1986])
at the University of California, San Francisco, is parameterized and defined only for proteins and DNA.
For the latter classes of molecules, various authors have added parameters and extended AMBER in other ways to suit their calculations.
The AMBER forcefield has also been made specifically applicable to polysaccharides (Homans 1990[^homans1990]).

The AMBER Force field has been extended with GAFF. In the following section you can see the most common atom types in AMBER and GAFF.
If your required atom type is not here, consult the documentation on GAFF (GAFF publication[^wang2004]).

For nucleic acids it is important to note that some atom types have been chosen due to the influence the aromatic rings have on their substituents.
Therefore some atoms have been assigned an orbital hybridisation type that is closer to their geometry than their actual hybridisation.
[Source: nonplanarity of NH2 groups](https://structbio.vanderbilt.edu/archives/amber-archive/2006/4057.php).
This information is important for when you are researching the practical side of a `.prep` file's atom typing.
Also the following publication is very important!
*Revisiting the planarity of nucleic acid bases: Pyramidilization at glycosidic nitrogen in purine bases is modulated by orientation of glycosidic torsion*[^sychrovsky2009]

### Atom types

#### Hydrogen types

*AMBER* types:

- **H**: amide or imino hydrogen
- **HC**: explicit hydrogen attached to aliphatic carbon
- **HO**: hydrogen on oxygen (hydroxyl groups)
- **HS**: hydrogen on sulfur
- **HW**: hydrogen in water
- **H2**: amino hydrogen in NH₂
- **H3**: hydrogen of lysine or arginine (positively charged)

*GAFF* types:

- **H1**: hydrogen on aliphatic carbon with 1 electron-withdrawal group
- **H2**: hydrogen on aliphatic carbon with 2 electron-withdrawal groups (this also coincides with H2)
- **H3**: hydrogen on aliphatic carbon with 3 electron-withdrawal groups (this also coincides with H3)
- **H4**: hydrogen on aromatic carbon with 1 electron-withdrawal group
- **H5**: hydrogen on aromatic carbon with 2 electron-withdrawal groups
- **HA**: explicit carbon attached to aromatic carbon
- **HN**: hydrogen on nitrogen (when the other hydrogen-nitrogen typing do not fit, generally)

#### All-atom carbon types

*AMBER* types:

- **C**: sp² carbonyl carbon and aromatic carbon with hydroxyl substituent in tyrosine
- **CA**: sp² aromatic carbon in 6-membered ring with 1 substituent
- **CB**: sp² aromatic carbon at junction between 5- and 6-membered rings
- **CC**: sp² aromatic carbon in 5-membered ring with 1 substituent and next to a nitrogen
- **CK**: sp² aromatic carbon in 5-membered ring between 2 nitrogens and bonded to 1 hydrogen (in purine)
- **CM**: sp² same as CJ but one substituent
- **CN**: sp² aromatic junction carbon in between 5- and 6-membered rings
- **CQ**: sp² carbon in 6-membered ring of purine between 2 NC nitrogens and bonded to 1 hydrogen
- **CR**: sp² aromatic carbon in 5-membered ring between 2 nitrogens and bonded to 1 H (in his)
- **CT**: sp³ carbon with 4 explicit substituents
- **CV**: sp² aromatic carbon in 5-membered ring bonded to 1 N and bonded to an explicit hydrogen
- **CW**: sp² aromatic carbon in 5-membered ring bonded to 1 N-H and bonded to an explicit hydrogen
- **C\***: sp² aromatic carbon in 5-membered ring with 1 substituent

#### United-atom carbon types

*AMBER* types:

- **CD**: sp² aromatic carbon in 6-membered ring with 1 hydrogen
- **CE**: sp² aromatic carbon in 5-membered ring between 2 nitrogens with 1 hydrogen (in purines)
- **CF**: sp² aromatic carbon in 5-membered ring next to a nitrogen without a hydrogen
- **CG**: sp² aromatic carbon in 5-membered ring next to an N-H
- **CH**: sp² carbon with 1 hydrogen
- **CI**: sp³ This atom type is also created since Pérez et al.,
  *Refinement of the AMBER Force Field for Nucleic Acids: Improving the description of α/γ conformers*[^perez2007], to assign to C5' carbons
- **CJ**: sp² carbon in pyrimidine at positions 5 or 6 (more pure double bond than aromatic with 1 hydrogen)
- **CP**: sp² aromatic carbon in 5-membered ring between 2 nitrogens with one hydrogen (in his)
- **C2**: sp² carbon with 2 hydrogens
- **C3**: sp² carbon with 3 hydrogens

*GAFF* types:

- **C3**: sp³ carbon (according to GAFF this can also be this)

#### Nitrogen types

*AMBER* types:

- **N**: sp² nitrogen in amide group
- **NA**: sp² nitrogen in 5-membered ring with hydrogen attached
- **NB**: sp² nitrogen in 5-membered ring with lone pairs
- **NC**: sp² nitrogen in 6-membered ring with lone pairs
- **NT**: sp² nitrogen with 3 substituents
- **N2**: sp² nitrogen in base NH₂ group or arginine NH₂
- **N3**: sp² nitrogen with 4 substituents
- **N\***: sp² nitrogen in purine or pyrimidine with alkyl group attached

*GAFF* types:

- **N3**: sp³ nitrogen with 3 substituents (according to GAFF this can also be this)
- **NH**: amine nitrogen connected to aromatic rings

#### Oxygen types

*AMBER* types:

- **O**: carbonyl oxygen
- **OH**: alcohol oxygen
- **OS**: ether or ester oxygen
- **OW**: water oxygen
- **O2**: carboxyl or phosphate nonbonded oxygen

#### Sulfur types

*AMBER* types:

- **S**: sulfur in disulfide linkage or methionine
- **SH**: sulfur in cystine

#### Phosphorus types

*AMBER* types:

- **P**: phosphorus in phosphate group

#### Ion types

*AMBER* types:

- **CU**: copper ion (Cu²⁺)
- **CO**: calcium ion (Ca²⁺)
- **I**: iodine ion (I⁻)
- **CL**: chlorine ion (Cl⁻)
- **MG**: magnesium ion (Mg²⁺)
- **QC**: cesium ion (Cs⁺)
- **QK**: potassium ion (K⁺)
- **QL**: lithium ion (Li⁺)
- **QN**: sodium ion (Na⁺)
- **QR**: rubidium ion (Rb⁺)

#### Others

*AMBER* types:

- **LP**: Lone pair

## RESP DocString

```text
source resp script : /path/to/amber/amber18/AmberTools/src/etc/resp.F

AMBER MAILING LIST MAILS CONCERNING TECHNICAL STUFF OF RESP
https://structbio.vanderbilt.edu/archives/amber-archive/2005/2429.php

https://structbio.vanderbilt.edu/archives/amber-archive/2005/2434.php

https://structbio.vanderbilt.edu/archives/amber-archive/2005/2440.php

some information concerning resp documentation : https://upjv.q4md-forcefieldtools.org/RED/resp/

    RESP   version 2.1     October 1994 Jim Caldwell
    RESP   version 2.0     September 1992
       Author: Christopher Bayly

    ESPFIT version 1.0 modified by Ian Gould to run in conjunction
           with gaussian 90.

    ESPFIT version 1.0 (G80UCSF):

              U.CHANDRA SINGH AND P.A.KOLLMAN

    All authors:

              DEPARTMENT OF PHARMACEUTICAL CHEMISTRY
              SCHOOL OF PHARMACY
              UNIVERSITY OF CALIFORNIA
              SAN FRANCISCO   CA 94143

    This current file has been specifically editted by Jérôme Rihon
    All information have been retrieved from the resp.F file



-------------------------------------------------------------------------------------------------------
|     CAVEAT : Notice all formatting strings and leading or trailing blank lines/spaces               |
|             THIS IS SUPER DUPER IMPORTANT! DO NOT SKIP THIS INFORMATION                             |
-------------------------------------------------------------------------------------------------------



  file        flag            fortran                     purpose
  name                        unit
  -----------------------------------------------------------------------------------------------------

  input       -i              5       required            input options

  output      -o              6       always produced     output of results

  punch       -p              7       always produced     synopsis of results

  qin         -q              3       optional            replacement charges
                                                          (restraints)

  qout        -t              19      always produced     output of current charges
                                                          (restraint that can be used later)

  espot       -e              10      required            input of MEP and coordinates
                                                          (note: these must be in atomic units,
                                                            not aengstrom)


  qwts        -w              4       optional            input of new weight factors

  esout       -s              20      optional            generated MEP values for new charges

 -----------------------------------------------------------------------------------------------------
```

```text
 RESP INPUT FILE SPECIFIC INFORMATION
 ---------------------------------------------------------------------

      THIS PROGRAM FITS THE QUANTUM MECHANICALLY CALCULATED
      POTENTIAL AT MOLECULAR SURFACES USING AN ATOM-CENTERED
      POINT CHARGE MODEL. THE MOLECULAR SURFACES ARE GENERATED
      BEYOND VANDER WAAL SURFACE IN ORDER TO MINIMISE OTHER
      CONTRIBUTIONS SUCH AS EXCHANGE REPULSION AND CHARGE TRANSFER

 ---------------------------------------------------------------------

      -1st-   TITLE       FORMAT(10A8)

 ---------------------------------------------------------------------

      -2nd-

            OPTIONS FOR THE JOB begin with " &cntrl"
                                end with   " &end"
            note leading blanks !!!!!!!!!!


        INOPT   =  0  ... NORMAL RUN
                =  1  ... CYCLE THROUGH A LIST OF DIFFERENT qwt
                                 read from -w unit
                                 This is an optional flag

        IOUTOPT =  0   NORMAL RUN
                =  1   write restart info of new esp etc to
                               unit -s (esout unit)

        IQOPT   =  0  ... use the q's which are read the -i unit
                =  1  ... RESET ALL INITIAL CHARGES TO ZERO
                =  2  ... READ IN NEW INITIAL CHARGES FROM -q (qwt)
                =  3  ... READ IN NEW INITIAL CHARGES FROM  -q (qwt)
                                   AND PERFORM AVERAGING OF THOSE NEW
                                   INITIAL CHARGES ACCORDING TO IVARY VALUES

        ihfree  =  0  ... ALL ATOMS ARE RESTRAINED
                =  1  ... HYDROGENS NOT RESTRAINED
                          Only heavy atoms are restrained

        irstrnt =  0  ... HARMONIC RESTRAINTS (old style)
                =  1  ... HYPERBOLIC RESTRAINT TO CHARGE OF ZERO (default)
                           Not recommend to do anything else than (1)
                =  2  ... ONLY ANALYSIS OF INPUT CHARGES; NO CHARGE FITTING IS CARRIED OUT

        qwt     =  restraint weight if irstrnt = 1 (default)
                =  0.0005 ... Recommend first stage restraints
                =  0.001  ... Recommend second stage restraints



       (DEPRECATED, IT IS ALWAYS IN ATOMIC UNITS [BOHR])
       (DO NOT MIND THIS VARIABLE)
       iunits    =  0  ... atom coordinates in angstroms
                 =  1       "     "          "  bohrs

  NOTE: ESP coordinates must always be in Bohrs

 --------------------------------------------------------------------------

      -3rd-
            wtmol .... relative weight for the molecule if multiple molecule fit (1.0 otherwise)

              FORMAT(F10.5) (10.5f - float - : ten characters of space allowed
                leading up to the decimal points, with 5 characters trailing)

 --------------------------------------------------------------------------

      -4th-
            subtitle for molecule

```

```text
 --------------------------------------------------------------------------

      -5th-
            CHARGE,IUNIQ ( THE NUMBER OF UNIQUE CENTERS for this molecule)
                         ( The charge of the total molecule , the amount of atoms in the molecule)

            FORMAT(2I5) (2 times 5d - int int - : five characters of space allowed for the integer value.)

 --------------------------------------------------------------------------
      -6th-

            Name, IVARY

        NAME    = ATOMIC number

        IVARY   = CONTROL OF CHARGE VARIATION OF THIS CENTER
                =  0 CHARGE VARIED INDEPENDENTLY OF PREVIOUS CENTERS
                         Let the charge of that atom (center) vary freely

                = -1 CHARGE FROZEN AT "INITIAL CHARGE" VALUE
                     Read the charge in from the file prompted to -q

                =  n CHARGE FITTED TOGETHER WITH CENTER n
                     Equivalence a charge with a similar atom (usually hydrogens in 2nd stage)
                                                              (hydrogens in amine in 1st stage)


            FORMAT(I5,i5) (2 times 5d - int int - : five characters of space allowed for
                the integer value.)
 -------------------------------------------------------------------------

      -7th-

            INTRA MOLECULE CHARGE CONSTRAINTS ...  blank line if no constraints

        ngrp    = number of centers in the group associated with this
                       constraint (i.e. the number of centers to be read in from `-q qin`)

        grpchg(i) = charge to which the associated group of atoms
                            (given on the next card) is to be contrained

            FORMAT(I5,F10.5) (5d 10.5f - int float - : five characters of space allowed
                for the integer value, followed by ten leading spaces then decimalpoint
                then five spaces for characters)

      -7.1-
        imol,iatom (16I5) (over multiple lines if there are more than eight (8) constraints.
                Not necessary to fill up the lines, only when needed)
                (16 times (per line) 5d - int x16 - : five characters of space allowed per integer value)

                the list (ngrp long) of the atom indices of those atoms to be
                constrained to the charge specified on the previous card.


    blank to end



 ------------------------------------------------------------------------
   -8th-

            INTERMOLECULAR CHARGE CONSTRAINTS  (atoms must sum to the specified value)

        same format as indvidual molecule constraints


    blank to end

 ------------------------------------------------------------------------

  -9th-

       MULTIPLE MOLECULE CONSTRAINTS .... constrain atoms on i to be the same as on j
       NGRP (I5) number of constaints
            (5d - int - : five characters of space allowed)
       (imol,iatom) (16I5) (over multiple lines if there are more than eight (8) constraints)
            (16 times (per line) 5d - int x16 - : five characters of space allowed per integer value )

       blank to end
```

```text
 THE FOLLOWING ARE SPECIFIC FORMATTING OF INPUTFILE
 --------------------------------------------------------------

    -e esp.dat

    Unit 10 input of ESP's  (mandatory)


    FILE ____________
        |    natoms,nesp (2i5)
        |            X , Y , Z  .   FORMAT (17X,3E16.7)
        |    QUPOT , X , Y , Z  .   FORMAT (1X,4E16.7)
        |____________


    QUPOT = THE QUANTUM MECHANICAL ELECTROSTATIC POTENTIAL ( A.U )

    X,Y,Z = THE COORDINATE AT WHICH THE POTENTIAL IS CALCULATED ( A.U )

    (DEPRECATED ) NOTE : THE PROGRAM G80UCSF WRITES IN THIS FORMAT
        BUT THE OUTPUT OF G90 MUST BE TRANSLATED (PROGRAM BOHR).

    (17X means 17 blank spaces. The reason is that these lines are the coordinates of the file
            and the first column will contain the espot values)
    (3E16.7 is expected three values per line to be formatted in a float-scientific notation.
            I did it in just float values and it worked fine.)
 ------------------------------------------------------------------------

   -q qin (optional)

    Unit 3 (qin) input of replacement charges if requested (this is done with the -q flag)
                        (over multiple lines if there are more than eight (8) constraints)

    A charge-input file that contains charges you want atoms to be restrained to

    Charges will only be read from the qin file if the following flag has been set to either value
    iqopt = 2
          = 3

        (8 f10.6)  (i = 1,iuniq)

 -------------------------------------------------------------

    -w weights (optional)

    Specific weights for the atoms in a molecule
    Unit 4 input if new weight factors if requested

        (i5)  nqwt  number of new weights to cycle thru
        (f10.5)  new weights (nqwt lines)


 COMMANDLINE PROMPT USAGE
 --------------------------------------------------------------


     usage: resp -i input -o output -p punch -q qin -t qout -e espot

                (optional) -w qwts -s esout



FORTRAN FORMAT (courtesy of https://upjv.q4md-forcefieldtools.org/RED/resp/)
 --------------------------------------------------------------
 2I5 format             2 blocks of 5 characters (I = integer)

 16I5 format           16 blocks of 5 characters (I = integer)

 F10.5 format           A block of 10 characters long (F = float) with 5 digits after the decimal point

 17X,3E16.7 format      A block of 17 space characters long followed by 3 blocks of 16 characters long
                        with 7 digits after the decimal point (using the scientific E notation)

 1X,4E16.7 format       A space character followed by 4 blocks of 16 characters long
                        with 7 digits after the decimal point (using the scientific E notation)
```

[^weiner1984]: Weiner SJ, Kollman PA, Case DA, Singh UC, Ghio C, Alagona G, Profeta S, Weiner P.
    A new force field for molecular mechanical simulation of nucleic acids and proteins.
    *J. Am. Chem. Soc.* 1984, 106, 765–784. DOI: [10.1021/ja00315a051](https://doi.org/10.1021/ja00315a051)

[^weiner1986]: Weiner SJ, Kollman PA, Nguyen DT, Case DA.
    An all atom force field for simulations of proteins and nucleic acids.
    *J. Comput. Chem.* 1986, 7, 230–252. DOI: [10.1002/jcc.540070216](https://doi.org/10.1002/jcc.540070216)

[^homans1990]: Homans SW.
    A molecular mechanical force field for the conformational analysis of oligosaccharides:
    comparison of theoretical and crystal structures of Manα1-3Manβ1-4GlcNAc.
    *Biochemistry* 1990, 29, 9110–9118. DOI: [10.1021/bi00491a003](https://doi.org/10.1021/bi00491a003)

[^wang2004]: Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA.
    Development and testing of a general amber force field.
    *J. Comput. Chem.* 2004, 25, 1157–1174. DOI: [10.1002/jcc.20035](https://doi.org/10.1002/jcc.20035)

[^sychrovsky2009]: Sychrovský V, Foldynová-Trantírková S, Špačková N, Robeyns K, Van Meervelt L, Blankenfeldt W, Vokáčová Z, Šponer J, Trantírek L.
    Revisiting the planarity of nucleic acid bases: Pyramidilization at glycosidic nitrogen in purine bases is modulated by orientation of glycosidic torsion.
    *Nucleic Acids Res.* 2009, 37, 7321–7331. DOI: [10.1093/nar/gkp783](https://doi.org/10.1093/nar/gkp783)

[^perez2007]: Pérez A, Marchán I, Svozil D, Šponer J, Cheatham TE, Laughton CA, Orozco M.
    Refinement of the AMBER Force Field for Nucleic Acids: Improving the Description of α/γ Conformers.
    *Biophys. J.* 2007, 92, 3817–3829. DOI: [10.1529/biophysj.106.097782](https://doi.org/10.1529/biophysj.106.097782)
