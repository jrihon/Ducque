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
TABLE_NUCLEOTIDES = {
"DA" : [ DH + "dna_adenosine_2endo.json", DH + "dna_phosphate.json"],
"DG" : [ DH + "dna_guanosine_2endo.json", DH + "dna_phosphate.json"],
"DC" : [ DH + "dna_cytidine_2endo.json", DH + "dna_phosphate.json"],
"DT" : [ DH + "dna_thymidine_2endo.json", DH + "dna_phosphate.json"],
"RA" : [ DH + "rna_adenosine_3endo.json", DH + "rna_phosphate.json"],
"RC" : [ DH + "rna_cytidine_3endo.json", DH + "rna_phosphate.json"],
"RG" : [ DH + "rna_guanosine_3endo.json", DH + "rna_phosphate.json"],
"RU" : [ DH + "rna_uracil_3endo.json", DH + "rna_phosphate.json"],
"DDA" : [ DH + "b-homodna_adenosine_1-4chair.json", DH + "dna_phosphate.json"],
"DDG" : [ DH + "b-homodna_guanosine_1-4chair.json", DH + "dna_phosphate.json"],
"DDC" : [ DH + "b-homodna_cytidine_1-4chair.json", DH + "dna_phosphate.json"],
"DDT" : [ DH + "b-homodna_thymidine_1-4chair.json", DH + "dna_phosphate.json"],
"HA" : [ DH + "hna_adenosine_1-4chair.json", DH + "dna_phosphate.json"],
"HG" : [ DH + "hna_guanosine_1-4chair.json", DH + "dna_phosphate.json"],
"HC" : [ DH + "hna_cytidine_1-4chair.json", DH + "dna_phosphate.json"],
"HT" : [ DH + "hna_thymidine_1-4chair.json", DH + "dna_phosphate.json"],
"XA" : [ DH + "xyna_adenosine_3endo.json", DH + "dna_phosphate.json"],
"XG" : [ DH + "xyna_guanosine_3endo.json", DH + "dna_phosphate.json"],
"XC" : [ DH + "xyna_cytidine_3endo.json", DH + "dna_phosphate.json"],
"XU" : [ DH + "xyna_uracil_3endo.json", DH + "dna_phosphate.json"],
"2MA" : [ DH + "2-ome-rna_adenosine_3endo.json", DH + "dna_phosphate.json"],
"2MG" : [ DH + "2-ome-rna_guanosine_3endo.json", DH + "dna_phosphate.json"],
"2MC" : [ DH + "2-ome-rna_cytidine_3endo.json", DH + "dna_phosphate.json"],
"2MU" : [ DH + "2-ome-rna_uracil_3endo.json", DH + "dna_phosphate.json"],
"CA" : [ DH + "cena_adenosine_3endo.json", DH + "rna_phosphate.json"],
"CG" : [ DH + "cena_guanosine_3endo.json", DH + "rna_phosphate.json"],
"CC" : [ DH + "cena_cytidine_3endo.json", DH + "rna_phosphate.json"],
"CT" : [ DH + "cena_thymidine_3endo.json", DH + "rna_phosphate.json"],
"DXA" : [DH + "dxyna_adenosine_3endo.json", DH + "dna_phosphate.json"],
"DXG" : [DH + "dxyna_guanosine_3endo.json", DH + "dna_phosphate.json"],
"DXC" : [DH + "dxyna_cytidine_3endo.json", DH + "dna_phosphate.json"],
"DXT" : [DH + "dxyna_thymidine_3endo.json", DH + "dna_phosphate.json"],
"MA" : [DH + "mna_adenosine_1-4chair.json", DH + "dna_phosphate.json"],
"MG" : [DH + "mna_guanosine_1-4chair.json", DH + "dna_phosphate.json"],
"MC" : [DH + "mna_cytidine_1-4chair.json", DH + "dna_phosphate.json"],
"MT" : [DH + "mna_thymidine_1-4chair.json", DH + "dna_phosphate.json"],
"TA" : [DH + "tna_adenosine_3_4_twist.json", DH + "dna_phosphate.json"],
"TC" : [DH + "tna_cytidine_3_4_twist.json", DH + "dna_phosphate.json"],
"TG" : [DH + "tna_guanosine_3_4_twist.json", DH + "dna_phosphate.json"],
"TT" : [DH + "tna_thymidine_3_4_twist.json", DH + "dna_phosphate.json"],
"PTA" : [DH + "photna_adenosine_3_4_twist.json", DH + "dna_phosphate.json"],
"PTC" : [DH + "photna_cytidine_3_4_twist.json", DH + "dna_phosphate.json"],
"PTG" : [DH + "photna_guanosine_3_4_twist.json", DH + "dna_phosphate.json"],
"PTT" : [DH + "photna_thymidine_3_4_twist.json", DH + "dna_phosphate.json"],
}



# -----------------------------------------
#       COMPLEMENTARY REPOSITORY
# -----------------------------------------
TABLE_CONFORMATIONS = {
"DA": [ DH + "dna_adenosine_2endo.json", DH + "dna_adenosine_3endo.json"],
"DC": [ DH + "dna_cytidine_2endo.json", DH + "dna_cytidine_3endo.json"],
"DG": [ DH + "dna_guanosine_2endo.json", DH + "dna_guanosine_3endo.json"],
"DT": [ DH + "dna_thymidine_2endo.json", DH + "dna_thymidine_3endo.json"],
"RA": [ DH + "rna_adenosine_3endo.json"],
"RC": [ DH + "rna_cytidine_3endo.json"],
"RG": [ DH + "rna_guanosine_3endo.json"],
"RU": [ DH + "rna_uracil_3endo.json"],
"DDA" : [ DH + "b-homodna_adenosine_1-4chair.json"], # "b-homodna_adenosine_4-1chair.json"],
"DDG" : [ DH + "b-homodna_guanosine_1-4chair.json"], #"b-homodna_guanosine_4-1chair.json"],
"DDC" : [ DH + "b-homodna_cytidine_1-4chair.json"], #"b-homodna_cytosine_4-1chair.json"],
"DDT" : [ DH + "b-homodna_thymidine_1-4chair.json"], # "b-homodna_thymidine_4-1chair.json",],
"HA" : [ DH + "hna_adenosine_1-4chair.json"], # "hna_adenosine_4-1chair.json"],
"HG" : [ DH + "hna_guanosine_1-4chair.json"], # "hna_guanosine_4-1chair.json"],
"HC" : [ DH + "hna_cytidine_1-4chair.json"], # "hna_cytidine_4-1chair.json"],
"HT" : [ DH + "hna_thymidine_1-4chair.json"], # "hna_thymidine_4-1chair.json"],
"XA" : [ DH + "xyna_adenosine_3endo.json"],
"XG" : [ DH + "xyna_guanosine_3endo.json"],
"XC" : [ DH + "xyna_cytidine_3endo.json"],
"XU" : [ DH + "xyna_uracil_3endo.json"],
"2MA" : [ DH + "2-ome-rna_adenosine_3endo.json"],
"2MG" : [ DH + "2-ome-rna_guanosine_3endo.json"],
"2MC" : [ DH + "2-ome-rna_cytidine_3endo.json"],
"2MU" : [ DH + "2-ome-rna_uracil_3endo.json"],
"CA": [ DH + "cena_adenosine_3endo.json"],
"CG": [ DH + "cena_guanosine_3endo.json"],
"CC": [ DH + "cena_cytidine_3endo.json"],
"CT": [ DH + "cena_thymidine_3endo.json"],
"DXA" : [DH + "dxyna_adenosine_3endo.json"],
"DXG" : [DH + "dxyna_guanosine_3endo.json"],
"DXC" : [DH + "dxyna_cytidine_3endo.json"],
"DXT" : [DH + "dxyna_thymidine_3endo.json"],
"MA" : [DH + "mna_adenosine_1-4chair.json"],
"MG" : [DH + "mna_guanosine_1-4chair.json"],
"MC" : [DH + "mna_cytidine_1-4chair.json"],
"MT" : [DH + "mna_thymidine_1-4chair.json"],
"TA" : [DH + "tna_adenosine_3_4_twist.json"],
"TC" : [DH + "tna_cytidine_3_4_twist.json"],
"TG" : [DH + "tna_guanosine_3_4_twist.json"],
"TT" : [DH + "tna_thymidine_3_4_twist.json"],
"PTA" : [DH + "photna_adenosine_3_4_twist.json"],
"PTC" : [DH + "photna_cytidine_3_4_twist.json"],
"PTG" : [DH + "photna_guanosine_3_4_twist.json"],
"PTT" : [DH + "photna_thymidine_3_4_twist.json"],
}



# -----------------------------------------
#           BACKBONE REPOSITORY
# -----------------------------------------
TABLE_BACKBONE = {
"PHOSPHATE" : ["P"],
"DNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"RNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"B-HOMODNA" : ["O4'", "C4'", "C5'", "C6'", "O6'"],
"HNA" : ["O4'", "C4'", "C5'", "C6'", "O6'"],
"XYNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"2-OME-RNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"CENA" : ["O3'", "C3'", "C4'", "C7'", "O7'"],
"DXYNA" : ["O3'", "C3'", "C4'", "C5'", "O5'"],
"MNA" : ["N3'", "C4'", "C5'", "C6'", "O6'"],
"TNA" : ["O2'", "C2'", "C3'", "O3'"],
"PHOTNA" : ["O2'", "C2'", "C3'", "O3'", "CP3'"],
}



# -----------------------------------------
#           LINKER REPOSITORY
# -----------------------------------------
TABLE_LINKER_BACKBONE = {
"PHOSPHATE" : ["P", "OP2", "OP1"],
}
