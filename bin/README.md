# Daedalus

<p align='center'> There are no mistakes, only happy little accidents. </p>
<br />

 **A Software for the purpose of building native and synthetic nucleic acid duplexes.**
    

#### All rights reserved to the Laboratory of Medicinal Chemistry @Rega Institute of Medical Research, Herestraat 49, 3000 Leuven, Belgium. Katholieke Universiteit Leuven (KUL).
    Software created and written by Doctorandus Rihon, Jerome

    Supervisor      prof. dr. Eveline Lescrinier
    Co-promotors    prof. dr. Vitor Bernardes Pinheiro
                    prof. dr. Matheus Froeyen

    Acknowledgement:
                    dr. Rinaldo Wander Montalvao for his guidance on the fundamentals of linear algebra.
                    dr. Charles-Alexandre Mattelaer for his guidance on Quantum Mechanics and without his experimental work, 
		    	Daedalus could have never been conceived.


Daedalus has three (3) main functions:
- --Daedalus : the nucleic acid builder
- --randomise : returns a randomised sequence to the user. This can then be read by --Daedalus
- --transmute : converts a given pdb file to the correct json format.
- --xyz_pdb : converts a given xyz file to the proper pdb format.

## Authors

[@jrihon](https://www.github.com/jrihon)

  
## Documentation

this

### Functions:
- #### Daedalus:

        Usage - $ python main.py --Daedalus INPUTFILE
        The inputfile is read in and the sequence is built accordingly.

        At any given time, there are two (2) flags in total that should be involved in the Daedalus - Nucleic Acid Builder.
        
        --sequence SEQUENCE
            Only valid input in the file just a string of the required nucleotides, separated by a comma and then a space.
                Example: dT, dC, dA, dA, dC, dG, dG, dT, dA

        --complement COMPLEMENT
            The complement flag denotes the structure of the complementary strand
            The following strings are valid inputs : homo (homoduplex)    DNA (DNA complementary)    RNA (RNA complementary) 
            A list of nucleotides is also a valid input, if one wants to specify the complementary strand. 
                Example: --complement homo
                Example: --complement dT, dA, dC, dC, dG, dT, dT, dG, dA

- #### Randomise:

        Usage - $ python main.py --randomise INPUTFILE
        The inputfile is read in and a sequence is outputted according to the prompted parameters.

        At any given time, there are two (2) flags in total that should be involved in the randomisation at all time.
            There are three in total, but two are mutually exclusive. (--length, --sequence) 

        --chemistry CHEMISTRY
            The chemistry that defines the prompted nucleotide or the chemistry that is involved with said chemistry
                Example: --chemistry DNA
                Example: --chemistry DNA, RNA, TNA
            
        --length LENGTH
            The length, in amount of nucleotides, of the sequence the user wants to generate. Required to be an integer.
                Example: --length 30

        --sequence SEQUENCE
            Provide the file with a sequence. Only the bases are required.
            The values should be prompted in the give order, separated by a comma and then a space:
                Example: --sequence A, C, T, G, G, A, A, T, C, A


        NB: Both a single (SINGLE) or multiple (LIST) chemistries can be prompted
            - randomise a sequence's bases with a given chemistry : 
                --chemistry SINGLE
                --length LENGTH

            - randomise the chemistries of a given sequence with a set of prompted chemistries : 
                --chemistry LIST
                --sequence SEQUENCE

            - randomise a sequence's bases with a set of prompted chemistries :
                --chemistry LIST
                --length LENGTH


- #### Transmute:

        Usage - $ python main.py --transmute INPUTFILE
        The inputfile is read in and the json file is formatted accordingly.

        There are five (5) flags total involved in the transmutation of a pdb structure to a json file
        --pdb PDB
            The name of the file of the structure you want to convert to json
                Example: --pdb dna_A.pdb

        --id ID
            The chemistry that defines the prompted nucleotide or the chemistry that is involved with said chemistry
                Example: --id DNA

        --comformation CONFORMATION
            The conformation that denotes the nucleic acid. Used to name the output .json file. Used when building the complementary strand.
            NB : for the leading strand, I typically use a conformation that is the most stable.
                The complementary strand is then fitted onto the leading strand.
                Example: --conformation 2endo
                Example: --conformation 1-4boat

        --moiety MOIETY
            The moiety that the structure defines. Should either be "nucleoside" or "linker"
                Example: --moiety nucleoside

        --dihedrals ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, CHI
            The dihedrals that involve the backbone and the anomeric carbon
            The values should be prompted in the give order, separated by a comma and then a space:
                Example: --dihedrals -39.246, -151.431, 30.929, 156.517, 159.171, -98.922, -99.315

        --bondangles ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, CHI
            The bond angles that involve the backbone and the anomeric carbon 
            The values should be prompted in the give order, separated by a comma and then a space:
                Example: --bondangles 101.407, 118.980, 110.017, 115.788, 111.943, 119.045, 126.013


- #### XYZ_PDB:
        Usage - $ python main.py --xyz_pdb INPUTFILE
        The inputfile is read, the xyz file used as an input to output a well formatted pdb.

        There are three(3) flags total involved in the conversion of a xyz coordinate file to a pdb structure file.
        --xyz XYZ
            The name of the file of the molecule you want to convert to pdb
               Example --xyz dna_2endo.xyz

        --atomID ATOMID
            The identifier for the molecule, typically named the 'Residue name' column. Right before the 'Chain' column. Typically a three-letter code, but can also be two or one.
            NB : Daedalus does not allow custom nucleic acid chemistries with an atomID unequal to three!
                Example: --atomID dXY
                    (in this example, dXA is equal to deoxy Xylose nucleic acid with an adenin base)

        --atomname_list ATOMNAME_LIST
            The ordered list of atoms that belong in the 'Atom name' column in a pdb file. Typically the third column, after the atom numbers.
            The order needs to so that it follows the order of the atoms from the xyz file. 
            Daedalus has a built-in method to check whether the order is correct by element, but the responsability is with the end-user to see everything is correct.
                Example: --atomname_list O5', C5', H5'1, H5'2, C4' ··· , O3' 


### JSON structure:

    { 
	pdb_properties : { Coordinates : [X, Y, Z] ,
                       Shape : (rows, columns) ,
                       Atoms : [] ,
                       Symbol : [] 
                      }

	identity: [ Chemistry, Abbreviation, Residue name, Base ]

	angles : { dihedrals : { Alpha : X,
                             Beta : X,
                             Gamma : X,
                             Delta : X,
                             Epsilon : X,
                             Zeta : X,
                             Chi : X
                             }
             }

             { bond_angles : { Alpha : X, 
                               Beta : X,
                               Gamma : X,
                               Delta : X,
                               Epsilon : X,
                               Chi : X
                               }
             }
    }

### Adding new chemistries:

BEFORE ADDING A NEW CHEMISTRY
In transmute_func_tools.py :
    Add to the nucleoside_dict the type of chemistry it is. The key should be in all caps.
    Add to the linker_dict the type of linker it corresponds with. The key should be in all caps.

    This only serves the purpose of identifying the json file when opening it as the user. This information is not used in the generation of nucleic acid duplexes.

BEFORE GENERATING A DUPLEX WITH THE NEW CHEMISTRY
In labyrinth_func_tools3.py :
    Add to the codex_acidum_nucleicum the most stable conformation of the chemistry you're adding to the library.
    Add to the complementary_codex all the conformations that you have at your disposal of the chemistry you're adding to the library.
        The key in both these abbreviated name of the nucleic acid chemistry.
    Add to the backbone_codex the sugar linker backbone of the chemistry you're adding to the library.

## Python environment

To run this project, you will need to add the libraries to your **python env**

`NumPy`, `SciPy`, `Pandas`
#
`$ pip install numpy` | `$ conda install -c numpy `

`$ pip install pandas` | `$ conda install -c pandas `

`$ pip install scipy` | `$ conda install -c scipy `
#
Daedalus also employs `sys`, `os` and `json`. These are built-in libraries, so no need to install these additionally.
The other python scripts are accompagnied when installing **Daedalus**.
  
## Acknowledgements
This README.md has been written with the help of [readme.so](https://readme.so)

