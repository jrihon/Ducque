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




1. TRANSMUTE NUCLEOBASE

[X] Rework transmute system!
    - [X] Make class for sugars
    - [X] Make class for linkers
    - [X] Make class for nucleobase

[X] Make Transmute module
    - [X] Make tinp file for transmute nucleobases
        - pdb file for atomic structure => passed by user
            -> atomic coordinates  > parsed from pdb
            -> atom names          > parsed from pdb
            -> size [x,3]          > parsed from pdb
        - abbr. name of res. (key)      => passed by user; like 5-Methyl Cytosine -> 5MC or something
        - atoms to rotate by            => passed by user ; N1, C2, C6
        - chemistry : name of the residue implemented. This becomes the filename.json 
    - [ ] GUI : Make --tbase flag to transmute nucleobases
    - [X] Allow --moiety `nucleobase` 

        "--pdb"
        "--chemistry"
        "--moiety"
        "--atoms"

    - [ ] For the Transmute GUI module : 
        - Write the queries out to a file 
        - Keep a counter on the GUI of how many modifications have been written to the file 
        - If by accident, the same query is written to the file, check against the values and disregard the write command if similar
            - When written to the file, clear all queries
        - Do not spawn extra elements in the GUI, make do with what is present, will work just fine





2. MODIFICATION
- make tables for querying if nucleobases exist
- parse pdb file to see if it some of the queries are viable
- iterate over pdb file make Pdb() instances of the nucleobases
- Copy all the pdb data and make into a class PdbFragment() or something
- Take note of the indices at which the residues that are going to be modified start and end (or their size) 

[X] Make table in ducquelib/library.p for parsing nucleobases => TABLE_NUCLEOBASE_MODS

[ ] Make PdbFragment() class to instance the specific residue that needs to be modified.
    - [X] Check if queries against PdbFile prompted is valid
    - [X] Only parse the residues in need of modification
    - [X] Apply modification only on those specific residues
    - [ ] Write out new file => DONT inplace modify the file

