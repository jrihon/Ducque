import sysDaedalus
""" This script will function as a repository for all the atoms, angles and dihedrals that need to be parsed to build up the nucleic acid duplex"""

# LEAVE THE `dh` VARIABLE ALONE
dh = sysDaedalus.return_DAEDALUS_home()
# LEAVE THE `dh` VARIABLE ALONE





# -----------------------------------------
#           NUCLEOSIDE REPOSITORY
# -----------------------------------------
codex_acidum_nucleicum = {
"dA" : [ dh + "json/dna_adenosine_2endo.json", dh + "json/dna_phosphate.json"],
"dG" : [ dh + "json/dna_guanosine_2endo.json", dh + "json/dna_phosphate.json"],
"dC" : [ dh + "json/dna_cytidine_2endo.json", dh + "json/dna_phosphate.json"],
"dT" : [ dh + "json/dna_thymidine_2endo.json", dh + "json/dna_phosphate.json"],
"rA" : [ dh + "json/rna_adenosine_3endo.json", dh + "json/rna_phosphate.json"],
"rC" : [ dh + "json/rna_cytidine_3endo.json", dh + "json/rna_phosphate.json"],
"rG" : [ dh + "json/rna_guanosine_3endo.json", dh + "json/rna_phosphate.json"],
"rU" : [ dh + "json/rna_uracil_3endo.json", dh + "json/rna_phosphate.json"],
"ddA" : [ dh + "json/b-homodna_adenosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"ddG" : [ dh + "json/b-homodna_guanosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"ddC" : [ dh + "json/b-homodna_cytidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"ddT" : [ dh + "json/b-homodna_thymidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hA" : [ dh + "json/hna_adenosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hG" : [ dh + "json/hna_guanosine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hC" : [ dh + "json/hna_cytidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"hT" : [ dh + "json/hna_thymidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"xA" : [ dh + "json/xylo_adenosine_3endo.json", dh + "json/dna_phosphate.json"],
"xG" : [ dh + "json/xylo_guanosine_3endo.json", dh + "json/dna_phosphate.json"],
"xC" : [ dh + "json/xylo_cytidine_3endo.json", dh + "json/dna_phosphate.json"],
"xU" : [ dh + "json/xylo_uracil_3endo.json", dh + "json/dna_phosphate.json"],
"2MA" : [ dh + "json/2-ome-rna_adenosine_3endo.json", dh + "json/dna_phosphate.json"],
"2MG" : [ dh + "json/2-ome-rna_guanosine_3endo.json", dh + "json/dna_phosphate.json"],
"2MC" : [ dh + "json/2-ome-rna_cytidine_3endo.json", dh + "json/dna_phosphate.json"],
"2MU" : [ dh + "json/2-ome-rna_uracil_3endo.json", dh + "json/dna_phosphate.json"],
"cA" : [ dh + "json/cena_adenosine_3endo.json", dh + "json/rna_phosphate.json"],
"cG" : [ dh + "json/cena_guanosine_3endo.json", dh + "json/rna_phosphate.json"],
"cC" : [ dh + "json/cena_cytidine_3endo.json", dh + "json/rna_phosphate.json"],
"cT" : [ dh + "json/cena_thymidine_3endo.json", dh + "json/rna_phosphate.json"],
"dxA" : [dh + "json/dxylo_adenosine_3endo.json", dh + "json/dna_phosphate.json"],
"dxG" : [dh + "json/dxylo_guanosine_3endo.json", dh + "json/dna_phosphate.json"],
"dxC" : [dh + "json/dxylo_cytidine_3endo.json", dh + "json/dna_phosphate.json"],
"dxT" : [dh + "json/dxylo_thymidine_3endo.json", dh + "json/dna_phosphate.json"],
#"mA" : [dh + "json/mna_adenosine_1-4chair.json", dh + "json/dna_phosphate.json"],
#"mG" : [dh + "json/mna_guanosine_1-4chair.json", dh + "json/dna_phosphate.json"],
#"mC" : [dh + "json/mna_cytidine_1-4chair.json", dh + "json/dna_phosphate.json"],
#"mT" : [dh + "json/mna_thymidine_1-4chair.json", dh + "json/dna_phosphate.json"],
"mA" : [dh + "json/mna_adenosine_2-4skew.json", dh + "json/dna_phosphate.json"],
"mG" : [dh + "json/mna_guanosine_2-4skew.json", dh + "json/dna_phosphate.json"],
"mC" : [dh + "json/mna_cytidine_2-4skew.json", dh + "json/dna_phosphate.json"],
"mT" : [dh + "json/mna_thymidine_2-4skew.json", dh + "json/dna_phosphate.json"],
}



# -----------------------------------------
#       COMPLEMENTARY REPOSITORY
# -----------------------------------------
conformations_codex = {
"dA": [ dh + "json/dna_adenosine_2endo.json", dh + "json/dna_adenosine_3endo.json"],
"dC": [ dh + "json/dna_cytidine_2endo.json", dh + "json/dna_cytidine_3endo.json"],
"dG": [ dh + "json/dna_guanosine_2endo.json", dh + "json/dna_guanosine_3endo.json"],
"dT": [ dh + "json/dna_thymidine_2endo.json", dh + "json/dna_thymidine_3endo.json"],
"rA": [ dh + "json/rna_adenosine_3endo.json"],
"rC": [ dh + "json/rna_cytidine_3endo.json"],
"rG": [ dh + "json/rna_guanosine_3endo.json"],
"rU": [ dh + "json/rna_uracil_3endo.json"],
"ddA" : [ dh + "json/b-homodna_adenosine_1-4chair.json"], # "json/b-homodna_adenosine_4-1chair.json"],
"ddG" : [ dh + "json/b-homodna_guanosine_1-4chair.json"], #"json/b-homodna_guanosine_4-1chair.json"],
"ddC" : [ dh + "json/b-homodna_cytidine_1-4chair.json"], #"json/b-homodna_cytosine_4-1chair.json"],
"ddT" : [ dh + "json/b-homodna_thymidine_1-4chair.json"], # "json/b-homodna_thymidine_4-1chair.json",],
"hA" : [ dh + "json/hna_adenosine_1-4chair.json"], # "json/hna_adenosine_4-1chair.json"],
"hG" : [ dh + "json/hna_guanosine_1-4chair.json"], # "json/hna_guanosine_4-1chair.json"],
"hC" : [ dh + "json/hna_cytidine_1-4chair.json"], # "json/hna_cytidine_4-1chair.json"],
"hT" : [ dh + "json/hna_thymidine_1-4chair.json"], # "json/hna_thymidine_4-1chair.json"],
"xA" : [ dh + "json/xylo_adenosine_3endo.json"],
"xG" : [ dh + "json/xylo_guanosine_3endo.json"],
"xC" : [ dh + "json/xylo_cytidine_3endo.json"],
"xU" : [ dh + "json/xylo_uracil_3endo.json"],
"2MA" : [ dh + "json/2-ome-rna_adenosine_3endo.json"],
"2MG" : [ dh + "json/2-ome-rna_guanosine_3endo.json"],
"2MC" : [ dh + "json/2-ome-rna_cytidine_3endo.json"],
"2MU" : [ dh + "json/2-ome-rna_uracil_3endo.json"],
"cA": [ dh + "json/cena_adenosine_3endo.json"],
"cG": [ dh + "json/cena_guanosine_3endo.json"],
"cC": [ dh + "json/cena_cytidine_3endo.json"],
"cT": [ dh + "json/cena_thymidine_3endo.json"],
"dxA" : [dh + "json/dxylo_adenosine_3endo.json"],
"dxG" : [dh + "json/dxylo_guanosine_3endo.json"],
"dxC" : [dh + "json/dxylo_cytidine_3endo.json"],
"dxT" : [dh + "json/dxylo_thymidine_3endo.json"],
#"mA" : [dh + "json/mna_adenosine_1-4chair.json"],
#"mG" : [dh + "json/mna_guanosine_1-4chair.json"],
#"mC" : [dh + "json/mna_cytidine_1-4chair.json"],
#"mT" : [dh + "json/mna_thymidine_1-4chair.json"],
"mA" : [dh + "json/mna_adenosine_2-4skew.json"],
"mG" : [dh + "json/mna_guanosine_2-4skew.json"],
"mC" : [dh + "json/mna_cytidine_2-4skew.json"],
"mT" : [dh + "json/mna_thymidine_2-4skew.json"],
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
