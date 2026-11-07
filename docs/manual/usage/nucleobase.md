# Nucleobase modification

The **Nucleobase** module allows to modify an existing structure with a simple query.
The nucleobases need to be recognised by ducque in order to be modifiable.

## Run the command

```bash
Ducque --nbase nucleobase-modificaton.inputfile
```

```text title="nucleobase-modificaton.inputfile
--nucleobase
-position 15.A  (required, where chain is optional)
-mod Psi        (required)
-resname pU     (optional, takes on resname of the current residue in the original pdb)
-reorient       (optional)
```

- position : the residue number and chain of the residue to change
- mod : the modification is needs to receive. Psi, but natural basepairs also work (for mismatching e.g.)
- resname : the residue name the full molecule should get
- reorient : to flip the regular orientation of the nucleobase (from watson-crick to hoogsteens)

## Run with more modification

You can run multiple modifications per ducque call.

```text title="nucleobase-modificaton-multiple.inputfile
--nucleobase
-position 15.A
-mod Psi
-resname pU
-reorient
--nucleobase
-position 3.B
-mod A
-resname rA
```

## gameplan contents

```text
$ Ducque --bases INPUTFILE.txt
```

```text
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

## TODO

1. TRANSMUTE NUCLEOBASE

[X] Rework transmute system!
    - [X] Make class for sugars
    - [X] Make class for linkers
    - [X] Make class for nucleobase

[X] Make Transmute module
    - [X] Make tinp file for transmute nucleobases
        - pdb file for atomic structure => passed by user
            -> atomic coordinates > parsed from pdb
            -> atom names > parsed from pdb
            -> size [x,3] > parsed from pdb
        - abbr. name of res. (key) => passed by user; like 5-Methyl Cytosine -> 5MC or something
        - atoms to rotate by            => passed by user ; N1, C2, C6
        - chemistry : name of the residue implemented. This becomes the filename.json
    - [ ] GUI : Make --tbase flag to transmute nucleobases
    - [X] Allow --moiety `nucleobase`

```text
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
```

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
