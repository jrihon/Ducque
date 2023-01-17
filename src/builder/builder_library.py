import systemsDucque

""" This script will function as a repository for all the atoms, angles and dihedrals that need to be parsed to build up the nucleic acid duplex"""

# LEAVE THE `DH` VARIABLE ALONE
DH = systemsDucque.return_DUCQUEHOME()




# -----------------------------------------
#           NUCLEOSIDE REPOSITORY
# -----------------------------------------

# -------------------
#                   NUCLEOSIDE              ,             LINKER
# -------------------
codex_acidum_nucleicum = {
"dA" : [ DH + "json/dna_adenosine_2endo.json", DH + "json/dna_phosphate.json"],
"dG" : [ DH + "json/dna_guanosine_2endo.json", DH + "json/dna_phosphate.json"],
"dC" : [ DH + "json/dna_cytidine_2endo.json", DH + "json/dna_phosphate.json"],
"dT" : [ DH + "json/dna_thymidine_2endo.json", DH + "json/dna_phosphate.json"],
"rA" : [ DH + "json/rna_adenosine_3endo.json", DH + "json/rna_phosphate.json"],
"rC" : [ DH + "json/rna_cytidine_3endo.json", DH + "json/rna_phosphate.json"],
"rG" : [ DH + "json/rna_guanosine_3endo.json", DH + "json/rna_phosphate.json"],
"rU" : [ DH + "json/rna_uracil_3endo.json", DH + "json/rna_phosphate.json"],
"ddA" : [ DH + "json/b-homodna_adenosine_1-4chair.json", DH + "json/dna_phosphate.json"],
"ddG" : [ DH + "json/b-homodna_guanosine_1-4chair.json", DH + "json/dna_phosphate.json"],
"ddC" : [ DH + "json/b-homodna_cytidine_1-4chair.json", DH + "json/dna_phosphate.json"],
"ddT" : [ DH + "json/b-homodna_thymidine_1-4chair.json", DH + "json/dna_phosphate.json"],
"hA" : [ DH + "json/hna_adenosine_1-4chair.json", DH + "json/dna_phosphate.json"],
"hG" : [ DH + "json/hna_guanosine_1-4chair.json", DH + "json/dna_phosphate.json"],
"hC" : [ DH + "json/hna_cytidine_1-4chair.json", DH + "json/dna_phosphate.json"],
"hT" : [ DH + "json/hna_thymidine_1-4chair.json", DH + "json/dna_phosphate.json"],
"xA" : [ DH + "json/xylo_adenosine_3endo.json", DH + "json/dna_phosphate.json"],
"xG" : [ DH + "json/xylo_guanosine_3endo.json", DH + "json/dna_phosphate.json"],
"xC" : [ DH + "json/xylo_cytidine_3endo.json", DH + "json/dna_phosphate.json"],
"xU" : [ DH + "json/xylo_uracil_3endo.json", DH + "json/dna_phosphate.json"],
"2MA" : [ DH + "json/2-ome-rna_adenosine_3endo.json", DH + "json/dna_phosphate.json"],
"2MG" : [ DH + "json/2-ome-rna_guanosine_3endo.json", DH + "json/dna_phosphate.json"],
"2MC" : [ DH + "json/2-ome-rna_cytidine_3endo.json", DH + "json/dna_phosphate.json"],
"2MU" : [ DH + "json/2-ome-rna_uracil_3endo.json", DH + "json/dna_phosphate.json"],
"cA" : [ DH + "json/cena_adenosine_3endo.json", DH + "json/rna_phosphate.json"],
"cG" : [ DH + "json/cena_guanosine_3endo.json", DH + "json/rna_phosphate.json"],
"cC" : [ DH + "json/cena_cytidine_3endo.json", DH + "json/rna_phosphate.json"],
"cT" : [ DH + "json/cena_thymidine_3endo.json", DH + "json/rna_phosphate.json"],
"dxA" : [DH + "json/dxylo_adenosine_3endo.json", DH + "json/dna_phosphate.json"],
"dxG" : [DH + "json/dxylo_guanosine_3endo.json", DH + "json/dna_phosphate.json"],
"dxC" : [DH + "json/dxylo_cytidine_3endo.json", DH + "json/dna_phosphate.json"],
"dxT" : [DH + "json/dxylo_thymidine_3endo.json", DH + "json/dna_phosphate.json"],
#"mA" : [DH + "json/mna_adenosine_1-4chair.json", DH + "json/dna_phosphate.json"],
#"mG" : [DH + "json/mna_guanosine_1-4chair.json", DH + "json/dna_phosphate.json"],
#"mC" : [DH + "json/mna_cytidine_1-4chair.json", DH + "json/dna_phosphate.json"],
#"mT" : [DH + "json/mna_thymidine_1-4chair.json", DH + "json/dna_phosphate.json"],
"mA" : [DH + "json/mna_adenosine_2-4skew.json", DH + "json/dna_phosphate.json"],
"mG" : [DH + "json/mna_guanosine_2-4skew.json", DH + "json/dna_phosphate.json"],
"mC" : [DH + "json/mna_cytidine_2-4skew.json", DH + "json/dna_phosphate.json"],
"mT" : [DH + "json/mna_thymidine_2-4skew.json", DH + "json/dna_phosphate.json"],
}



# -----------------------------------------
#       COMPLEMENTARY REPOSITORY
# -----------------------------------------
conformations_codex = {
"dA": [ DH + "json/dna_adenosine_2endo.json", DH + "json/dna_adenosine_3endo.json"],
"dC": [ DH + "json/dna_cytidine_2endo.json", DH + "json/dna_cytidine_3endo.json"],
"dG": [ DH + "json/dna_guanosine_2endo.json", DH + "json/dna_guanosine_3endo.json"],
"dT": [ DH + "json/dna_thymidine_2endo.json", DH + "json/dna_thymidine_3endo.json"],
"rA": [ DH + "json/rna_adenosine_3endo.json"],
"rC": [ DH + "json/rna_cytidine_3endo.json"],
"rG": [ DH + "json/rna_guanosine_3endo.json"],
"rU": [ DH + "json/rna_uracil_3endo.json"],
"ddA" : [ DH + "json/b-homodna_adenosine_1-4chair.json"], # "json/b-homodna_adenosine_4-1chair.json"],
"ddG" : [ DH + "json/b-homodna_guanosine_1-4chair.json"], #"json/b-homodna_guanosine_4-1chair.json"],
"ddC" : [ DH + "json/b-homodna_cytidine_1-4chair.json"], #"json/b-homodna_cytosine_4-1chair.json"],
"ddT" : [ DH + "json/b-homodna_thymidine_1-4chair.json"], # "json/b-homodna_thymidine_4-1chair.json",],
"hA" : [ DH + "json/hna_adenosine_1-4chair.json"], # "json/hna_adenosine_4-1chair.json"],
"hG" : [ DH + "json/hna_guanosine_1-4chair.json"], # "json/hna_guanosine_4-1chair.json"],
"hC" : [ DH + "json/hna_cytidine_1-4chair.json"], # "json/hna_cytidine_4-1chair.json"],
"hT" : [ DH + "json/hna_thymidine_1-4chair.json"], # "json/hna_thymidine_4-1chair.json"],
"xA" : [ DH + "json/xylo_adenosine_3endo.json"],
"xG" : [ DH + "json/xylo_guanosine_3endo.json"],
"xC" : [ DH + "json/xylo_cytidine_3endo.json"],
"xU" : [ DH + "json/xylo_uracil_3endo.json"],
"2MA" : [ DH + "json/2-ome-rna_adenosine_3endo.json"],
"2MG" : [ DH + "json/2-ome-rna_guanosine_3endo.json"],
"2MC" : [ DH + "json/2-ome-rna_cytidine_3endo.json"],
"2MU" : [ DH + "json/2-ome-rna_uracil_3endo.json"],
"cA": [ DH + "json/cena_adenosine_3endo.json"],
"cG": [ DH + "json/cena_guanosine_3endo.json"],
"cC": [ DH + "json/cena_cytidine_3endo.json"],
"cT": [ DH + "json/cena_thymidine_3endo.json"],
"dxA" : [DH + "json/dxylo_adenosine_3endo.json"],
"dxG" : [DH + "json/dxylo_guanosine_3endo.json"],
"dxC" : [DH + "json/dxylo_cytidine_3endo.json"],
"dxT" : [DH + "json/dxylo_thymidine_3endo.json"],
#"mA" : [DH + "json/mna_adenosine_1-4chair.json"],
#"mG" : [DH + "json/mna_guanosine_1-4chair.json"],
#"mC" : [DH + "json/mna_cytidine_1-4chair.json"],
#"mT" : [DH + "json/mna_thymidine_1-4chair.json"],
"mA" : [DH + "json/mna_adenosine_2-4skew.json"],
"mG" : [DH + "json/mna_guanosine_2-4skew.json"],
"mC" : [DH + "json/mna_cytidine_2-4skew.json"],
"mT" : [DH + "json/mna_thymidine_2-4skew.json"],
}



# -----------------------------------------
#           BACKBONE REPOSITORY
# -----------------------------------------
backbone_codex = {
"Phosphate" : ["P"],
"DNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"RNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"b-homoDNA" : ["O4'", "C4'", "C5'", "C6'", "O6'"],
"HNA" : ["O4'", "C4'", "C5'", "C6'", "O6'"],
"Xylo" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"2-OMe-RNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"CeNA" : ["O3'", "C3'", "C4'", "C7'", "O7'"],
"dXylo" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"MNA" : ["N3'", "C4'", "C5'", "C6'", "O6'"],
}



# -----------------------------------------
#           LINKER REPOSITORY
# -----------------------------------------
linker_codex = {
"Phosphate" : ["P", "OP2", "OP1"],
}
