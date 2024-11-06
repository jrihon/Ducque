# Ducque
The model builder of choice to generate helical structures. Provides an interface for the user to supply the Ducque library with user-defined, custom chemistries so you can build virtual structures with never-seen-before molecules!

Ducque has four (4) main functions:
- `--build` : module to generate (X)NA helical models.
- `--transmute` : module to insert user-defined chemistries into Ducque's library.
- `--xyz_pdb` : module to convert a given xyz file to the pdb format.
- `--randomise` : module to randomise a sequence, used for `--build` .
- `--gui` opens the GUI to use the modules up above, instead of the CLI. The `tlinker` module, in the `gui` module, is a specific preset for implementing linker fragments.

## Model building documentation
A comprehensive manual of Ducque is found inside the `Ducque/docs/` directory !

## Forcefield
To know more about how to build a forcefield, definitely check the `Ducque/ff/` directory!

### Charge derivation implementation with ORCA
All scripts to calculate point charges, according to the AMBER-compatible MK population analysis scheme, can be found in the `Ducque/ff/scripts/` directory.

The `How_To_Create_A_Forcefield.pdf` file guides the user through the correct protocol and where to download the requires software packages.
  
## Reference
If you've used Ducque or any of its parts, or you've used the `How_To_Create_A_Forcefield.pdf` as a guide for your own research, please cite the following published article : 

[1] Jérôme Rihon, Charles-Alexandre Mattelaer, Rinaldo Wander Montalvão, Mathy Froeyen, Vitor Bernardes Pinheiro, Eveline Lescrinier, Structural insights into the morpholino nucleic acid/RNA duplex using the new XNA builder Ducque in a molecular modeling pipeline, Nucleic Acids Research, 2024; DOI: 10.1093/nar/gkae135, [https://doi.org/10.1093/nar/gkae135](https://doi.org/10.1093/nar/gkae135)

```bibtex
@article{Rihon_2024,
    title={Structural insights into the morpholino nucleic acid/RNA duplex using the new XNA builder Ducque in a molecular modeling pipeline},
    volume={52},
    ISSN={1362-4962},
    url={http://dx.doi.org/10.1093/nar/gkae135},
    DOI={10.1093/nar/gkae135},
    number={6},
    journal={Nucleic Acids Research},
    publisher={Oxford University Press (OUP)},
    author={Rihon, Jérôme and Mattelaer, Charles-Alexandre and Montalvão, Rinaldo Wander and Froeyen, Mathy and Pinheiro, Vitor Bernardes and Lescrinier, Eveline},
    year={2024},
    month=feb,
    pages={2836–2847}
}
```
## Authors
Jérôme Rihon ([@jrihon](https://www.github.com/jrihon) )
