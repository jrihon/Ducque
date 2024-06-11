# inputs


```
$ Ducque --bases INPUTFILE.txt
```


```
--pdb FILE_WITH_STRUCTURES.pdb


--nucleobase 15.A , Psi, pU, HG

--nucleobase 
-position 15.A  (required, where chain is optional)
-mod Psi        (required)
-resname pU     (required)
-orientation HG (optional)

--nucleobase 
-position 2
-mod U
-resname U              => reasoning is that the software should not guess if the U modification comes from Uridine or the user specifies their own U modification
-orientation HG (optional)
```


# TODO

1. CLI PARSING
- make tables for querying if nucleobases exist
- parse pdb file to see if it some of the queries are viable


2. MODIFICATION
- iterate over pdb file make Pdb() instances of the nucleobases
- Copy all the pdb data and make into a class PdbFragment() or something
- Take note of the indices at which the residues that are going to be modified start and end (or their size) 
