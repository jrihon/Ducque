import systemsDucque

""" This script will function as a repository for all the atoms, angles and dihedrals that need to be parsed to build up the nucleic acid duplex"""

# LEAVE THE `DH` VARIABLE ALONE ; path/to/Ducque/json
DH = systemsDucque.return_DUCQUEHOME() + "json/"



# -----------------------------------------
#           NUCLEOSIDE REPOSITORY
# -----------------------------------------

# -------------------
#                   NUCLEOSIDE              ,             LINKER
# -------------------
codex_acidum_nucleicum = {
"dA" : [ DH + "dna_adenosine_2endo.json", DH + "dna_phosphate.json"],
"dG" : [ DH + "dna_guanosine_2endo.json", DH + "dna_phosphate.json"],
"dC" : [ DH + "dna_cytidine_2endo.json", DH + "dna_phosphate.json"],
"dT" : [ DH + "dna_thymidine_2endo.json", DH + "dna_phosphate.json"],
"rA" : [ DH + "rna_adenosine_3endo.json", DH + "rna_phosphate.json"],
"rC" : [ DH + "rna_cytidine_3endo.json", DH + "rna_phosphate.json"],
"rG" : [ DH + "rna_guanosine_3endo.json", DH + "rna_phosphate.json"],
"rU" : [ DH + "rna_uracil_3endo.json", DH + "rna_phosphate.json"],
"ddA" : [ DH + "b-homodna_adenosine_1-4chair.json", DH + "dna_phosphate.json"],
"ddG" : [ DH + "b-homodna_guanosine_1-4chair.json", DH + "dna_phosphate.json"],
"ddC" : [ DH + "b-homodna_cytidine_1-4chair.json", DH + "dna_phosphate.json"],
"ddT" : [ DH + "b-homodna_thymidine_1-4chair.json", DH + "dna_phosphate.json"],
"hA" : [ DH + "hna_adenosine_1-4chair.json", DH + "dna_phosphate.json"],
"hG" : [ DH + "hna_guanosine_1-4chair.json", DH + "dna_phosphate.json"],
"hC" : [ DH + "hna_cytidine_1-4chair.json", DH + "dna_phosphate.json"],
"hT" : [ DH + "hna_thymidine_1-4chair.json", DH + "dna_phosphate.json"],
"xA" : [ DH + "xylo_adenosine_3endo.json", DH + "dna_phosphate.json"],
"xG" : [ DH + "xylo_guanosine_3endo.json", DH + "dna_phosphate.json"],
"xC" : [ DH + "xylo_cytidine_3endo.json", DH + "dna_phosphate.json"],
"xU" : [ DH + "xylo_uracil_3endo.json", DH + "dna_phosphate.json"],
"2MA" : [ DH + "2-ome-rna_adenosine_3endo.json", DH + "dna_phosphate.json"],
"2MG" : [ DH + "2-ome-rna_guanosine_3endo.json", DH + "dna_phosphate.json"],
"2MC" : [ DH + "2-ome-rna_cytidine_3endo.json", DH + "dna_phosphate.json"],
"2MU" : [ DH + "2-ome-rna_uracil_3endo.json", DH + "dna_phosphate.json"],
"cA" : [ DH + "cena_adenosine_3endo.json", DH + "rna_phosphate.json"],
"cG" : [ DH + "cena_guanosine_3endo.json", DH + "rna_phosphate.json"],
"cC" : [ DH + "cena_cytidine_3endo.json", DH + "rna_phosphate.json"],
"cT" : [ DH + "cena_thymidine_3endo.json", DH + "rna_phosphate.json"],
"dxA" : [DH + "dxylo_adenosine_3endo.json", DH + "dna_phosphate.json"],
"dxG" : [DH + "dxylo_guanosine_3endo.json", DH + "dna_phosphate.json"],
"dxC" : [DH + "dxylo_cytidine_3endo.json", DH + "dna_phosphate.json"],
"dxT" : [DH + "dxylo_thymidine_3endo.json", DH + "dna_phosphate.json"],
"mA" : [DH + "mna_adenosine_1-4chair.json", DH + "dna_phosphate.json"],
"mG" : [DH + "mna_guanosine_1-4chair.json", DH + "dna_phosphate.json"],
"mC" : [DH + "mna_cytidine_1-4chair.json", DH + "dna_phosphate.json"],
"mT" : [DH + "mna_thymidine_1-4chair.json", DH + "dna_phosphate.json"],
#"mA" : [DH + "mna_adenosine_2-4skew.json", DH + "dna_phosphate.json"],
#"mG" : [DH + "mna_guanosine_2-4skew.json", DH + "dna_phosphate.json"],
#"mC" : [DH + "mna_cytidine_2-4skew.json", DH + "dna_phosphate.json"],
#"mT" : [DH + "mna_thymidine_2-4skew.json", DH + "dna_phosphate.json"],
"tA" : [DH + "tna_adenosine_3_4_twist.json", DH + "dna_phosphate.json"],
"tC" : [DH + "tna_cytidine_3_4_twist.json", DH + "dna_phosphate.json"],
"tG" : [DH + "tna_guanosine_3_4_twist.json", DH + "dna_phosphate.json"],
"tT" : [DH + "tna_thymidine_3_4_twist.json", DH + "dna_phosphate.json"],
"ptA" : [DH + "photna_adenosine_3_4_twist.json", DH + "dna_phosphate.json"],
"ptC" : [DH + "photna_cytidine_3_4_twist.json", DH + "dna_phosphate.json"],
"ptG" : [DH + "photna_guanosine_3_4_twist.json", DH + "dna_phosphate.json"],
"ptT" : [DH + "photna_thymidine_3_4_twist.json", DH + "dna_phosphate.json"],
}



# -----------------------------------------
#       COMPLEMENTARY REPOSITORY
# -----------------------------------------
conformations_codex = {
"dA": [ DH + "dna_adenosine_2endo.json", DH + "dna_adenosine_3endo.json"],
"dC": [ DH + "dna_cytidine_2endo.json", DH + "dna_cytidine_3endo.json"],
"dG": [ DH + "dna_guanosine_2endo.json", DH + "dna_guanosine_3endo.json"],
"dT": [ DH + "dna_thymidine_2endo.json", DH + "dna_thymidine_3endo.json"],
"rA": [ DH + "rna_adenosine_3endo.json"],
"rC": [ DH + "rna_cytidine_3endo.json"],
"rG": [ DH + "rna_guanosine_3endo.json"],
"rU": [ DH + "rna_uracil_3endo.json"],
"ddA" : [ DH + "b-homodna_adenosine_1-4chair.json"], # "b-homodna_adenosine_4-1chair.json"],
"ddG" : [ DH + "b-homodna_guanosine_1-4chair.json"], #"b-homodna_guanosine_4-1chair.json"],
"ddC" : [ DH + "b-homodna_cytidine_1-4chair.json"], #"b-homodna_cytosine_4-1chair.json"],
"ddT" : [ DH + "b-homodna_thymidine_1-4chair.json"], # "b-homodna_thymidine_4-1chair.json",],
"hA" : [ DH + "hna_adenosine_1-4chair.json"], # "hna_adenosine_4-1chair.json"],
"hG" : [ DH + "hna_guanosine_1-4chair.json"], # "hna_guanosine_4-1chair.json"],
"hC" : [ DH + "hna_cytidine_1-4chair.json"], # "hna_cytidine_4-1chair.json"],
"hT" : [ DH + "hna_thymidine_1-4chair.json"], # "hna_thymidine_4-1chair.json"],
"xA" : [ DH + "xylo_adenosine_3endo.json"],
"xG" : [ DH + "xylo_guanosine_3endo.json"],
"xC" : [ DH + "xylo_cytidine_3endo.json"],
"xU" : [ DH + "xylo_uracil_3endo.json"],
"2MA" : [ DH + "2-ome-rna_adenosine_3endo.json"],
"2MG" : [ DH + "2-ome-rna_guanosine_3endo.json"],
"2MC" : [ DH + "2-ome-rna_cytidine_3endo.json"],
"2MU" : [ DH + "2-ome-rna_uracil_3endo.json"],
"cA": [ DH + "cena_adenosine_3endo.json"],
"cG": [ DH + "cena_guanosine_3endo.json"],
"cC": [ DH + "cena_cytidine_3endo.json"],
"cT": [ DH + "cena_thymidine_3endo.json"],
"dxA" : [DH + "dxylo_adenosine_3endo.json"],
"dxG" : [DH + "dxylo_guanosine_3endo.json"],
"dxC" : [DH + "dxylo_cytidine_3endo.json"],
"dxT" : [DH + "dxylo_thymidine_3endo.json"],
"mA" : [DH + "mna_adenosine_1-4chair.json"],
"mG" : [DH + "mna_guanosine_1-4chair.json"],
"mC" : [DH + "mna_cytidine_1-4chair.json"],
"mT" : [DH + "mna_thymidine_1-4chair.json"],
#"mA" : [DH + "mna_adenosine_2-4skew.json"],
#"mG" : [DH + "mna_guanosine_2-4skew.json"],
#"mC" : [DH + "mna_cytidine_2-4skew.json"],
#"mT" : [DH + "mna_thymidine_2-4skew.json"],
"tA" : [DH + "tna_adenosine_3_4_twist.json"],
"tC" : [DH + "tna_cytidine_3_4_twist.json"],
"tG" : [DH + "tna_guanosine_3_4_twist.json"],
"tT" : [DH + "tna_thymidine_3_4_twist.json"],
"ptA" : [DH + "photna_adenosine_3_4_twist.json"],
"ptC" : [DH + "photna_cytidine_3_4_twist.json"],
"ptG" : [DH + "photna_guanosine_3_4_twist.json"],
"ptT" : [DH + "photna_thymidine_3_4_twist.json"],
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
"TNA" : ["O2'", "C2'", "C3'", "O3'"],
"PhoTNA" : ["O2'", "C2'", "C3'", "O3'", "CP3'"],
}



# -----------------------------------------
#           LINKER REPOSITORY
# -----------------------------------------
linker_codex = {
"Phosphate" : ["P", "OP2", "OP1"],
}
