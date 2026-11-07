# Installing and running Ducque

## Installation

The installation is done by cloning the repository, through the use of `git`.
It is advised to clone the repository to the desired location at the `$HOME`-level (`$HOME`, `$HOME/Desktop`, …).

Enter the Ducque directory and run the `setup.sh` script to give the correct permissions to the Ducque executable, located in `Ducque/bin`.
Finally, add the last line to your bash environment, by adding it to your `$HOME/.bashrc`.
MacOS users may replace mentions of `$HOME/.bashrc` with `$HOME/.bash_profile`.

For future reference, the `path/to/program` is the location on your system where Ducque is installed.
E.g., if Ducque is installed at `$HOME/Desktop/Ducque`, then `$HOME/Desktop` = `path/to/program`.

```bash
# Install on your system
git clone https://www.github.com/jrihon/Ducque.git # clone repo locally
cd Ducque

# Run the setup.sh script
chmod a+x setup.sh    # Give permission to run as an executable script
./setup.sh

# Add to ~/.bashrc
export PATH="$PATH:path/to/program/Ducque/bin"
```

## Environment

To install the correct packages for Ducque, it is advised to use the `conda` package manager for Python modules.
[Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

To use Ducque, it is also advised to create a separate virtual environment so as to not override or break any dependencies in the base environment.
To leave the virtual environment, the user can `deactivate` the conda environment to leave the **Ducque environment**.

The user must `activate` the Ducque environment before installing the dependencies.
The packages to install are **NumPy**, **SciPy** and **Tkinter**.

```bash
# Create a virtual environment using conda
conda create --name Ducque
# Enter the Ducque `conda environment`
conda activate Ducque

# Using the conda package manager install : NumPy, SciPy, Tkinter
conda install -c numpy scipy tk

# To leave the virtual environment
conda deactivate
```
