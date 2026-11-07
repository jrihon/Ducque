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
