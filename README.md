# Ducque

## Documentation
Ducque has four (4) main functions:
- `--build` : the nucleic acid builder
- `--transmute` : converts a given pdb file to the correct json format.
- `--xyz_pdb` : converts a given xyz file to the proper pdb format.
- `--randomise` : returns a randomised sequence to the user.
- `--gui` : open the GUI to use the modules up above, instead of the CLI

A comprehensive manual of Ducque is find inside the `Ducque/docs/` directory !

To know more about how to build a forcefield, definitely check out the `pdf` inside the `Ducque/ff/` directory!

<!--
## Environment

### Shell
To access the Ducque software from anywhere on your machine, add the following line to your `~/.bashrc` .<br />
Where path/to/program is the path to where you've installed Ducque. (`$ pwd`) inside the Ducque directory if you're unsure.
`export PATH=$PATH:path/to/program/Ducque/bin`

### Python
To run this project, you will need to add the libraries to your **python env** : `NumPy`, `SciPy`
#### Installation (conda or pip)
`$ pip install numpy` | `$ conda install -c numpy `

`$ pip install scipy` | `$ conda install -c scipy `
</br>
</br>
Ducque also employs `sys`, `os`, `tkinter` and `json`. These are built-in libraries, so no need to install these additionally.</br>
</br>
Depending on your `python3` version, you might have to install `tkinter` separately.
```shell
# Check if you have tkinter installed
$ python3 -m tkinter
```
-->
  
## Authors
Jérôme Rihon ([@jrihon](https://www.github.com/jrihon) )
