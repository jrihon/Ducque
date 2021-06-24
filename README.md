#### All rights reserved to the Laboratory of Medicinal Chemistry @Rega Institute of Medical Research, Herestraat 49, 3000 Leuven, Belgium. Katholieke Universiteit Leuven (KUL).
#    Software created and written by Doctorandus Rihon, Jerome
#
#    Supervisor      prof. dr. Eveline Lescrinier
#    Co-promotors    prof. dr. Vitor Bernardes Pinheiro
#                    prof. dr. Matheus Froeyen
#
#    Acknowledgement:
#                    dr. Rinaldo Wander Montalvao for his guidance on the fundamentals of linear algebra.
#                    dr. Charles-Alexandre Mattelaer for his guidance on Quantum Mechanics and without his experimental work, Daedalus could have never been conceived.
#
#
#
#                                                        " Nobody told me it wasn't feasible, so I went ahead and did it "
#
#
####

Daedalus has three (3) main functions:
    --Daedalus : the nucleic acid builder
    --randomise : returns a randomised sequence to the user. This can then be read by --Daedalus
    --transmute : converts a given pdb file to the correct json format.




#########
Daedalus:

    Usage - $ python main.py --Daedalus INPUTFILE
    The inputfile is read in and the sequence is built accordingly.

    Only valid inputs in the file just a string of the required nucleotides, separated by a comma and then a space :
        Example: dT, dC, dA, dA, dC, dG, dG, dT, dA




##########
Randomise:

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


##########
Transmute:

    Usage - $ python main.py --transmute INPUTFILE
    The inputfile is read in and the json file is formatted accordingly.

    There are five (5) flags total involved in the transmutation of a pdb structure to a json file
    --pdb PDB
        The name of the file of the structure you want to convert to json
            Example: --pdb dna_A.pdb

    --id ID
        The chemistry that defines the prompted nucleotide or the chemistry that is involved with said chemistry
            Example: --id DNA

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




##############
JSON structure:

{ 
	pdb_properties : { Coordinates : [X, Y, Z] ,
                       Shape : () ,
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
