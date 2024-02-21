# Ducque

Ducque has four (4) main functions:
- `--build` : the nucleic acid builder
- `--transmute` : converts a given pdb file to the correct json format.
- `--xyz_pdb` : converts a given xyz file to the proper pdb format.
- `--randomise` : returns a randomised sequence to the user.
- `--gui` opens the GUI to use the modules up above, instead of the CLI. The `tlinker` module, in the `gui` module, is a specific preset for implementing linker fragments.

## Documentation
A comprehensive manual of Ducque is found inside the `Ducque/docs/` directory !

## Forcefield
To know more about how to build a forcefield, definitely check the `Ducque/ff/` directory!

### Charge derivation implementation with ORCA
All scripts to calculate point charges, according to the AMBER-compatible MK population analysis scheme, can be found in the `Ducque/ff/scripts/` directory.
<br>
The `How_To_Create_A_Forcefield.pdf` file guides the user through the correct protocol and where to download the requires software packages.
  
## Authors
Jérôme Rihon ([@jrihon](https://www.github.com/jrihon) )
