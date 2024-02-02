import numpy as np
import os
import sys
import subprocess

"""

(c) 2022 Jérôme Rihon

License: MIT

Description: Generate a Connolly surface (Solvent Excluded Surface) around the molecule of interest and extract surface vertrices.
            Map the densities from the orbital with `orca_vpot` (Molecular Electrostatic Potential) onto the gridpoints.

            Congratulations, you have employed the Merz-Kollman population analysis scheme to generate ESP charges.

            This script then processes everything to conceive `resp` ready files.

            This script assumes you have installed the following binary (Linux only) : https://ccsb.scripps.edu/msms/downloads/
            And that you have made a copy of one of the two provided binaries and have copied it to a `msms` binary.




------------------------------------------------------------------------
Usage : $ python MK_PopAnalScheme.py template_Orca
    Where `template_Orca` is the name of the (HF - 6-31G*) outputfile generated.

------------------------------------------------------------------------
Arguments : As the first variable, give the basename of the file that was outputted when computing the HF 6-31G* single point energy.

        This script will search for `reoriented_*.pdb` in the current working directory to parse the atom names to assign vdW radii
        So make sure to have used the `reorient_nucleosides.py` before using this script!


"""

def parse_pdb_for_scraps() -> list:
    """ Atom name :       line 13 - 16 """

    pdb_prefix = "reoriented"

    cwd_list = [ i for i in os.listdir(os.getcwd()) if i.startswith(pdb_prefix) ]

    # if cwd_list is empty, stop the software
    if len(cwd_list) == 0:
        sys.exit("There are no files in the current working directory that has 'model' in the name of a pdb file.\n"
                "Make sure you have reoriented your molecules before parametrisation!\n")

    if isinstance(cwd_list, list):
        pdb_fname = cwd_list[0]
    else:
        pdb_fname = cwd_list


    # Check if file is in cwd or in the pdb directory
    try:
        os.path.isfile(pdb_fname)
    except FileNotFoundError:
        print(f"Could not find {pdb_fname} in the directory.\n")

    # Start new lists to append it all
    atomname = []

    # Read the file and fill out the dataframe
    with open(pdb_fname) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM' or line[:6] == 'HETATM':
                atomname.append(line[12:16])


    # Add the atom name list as an attribute
    atomname_list = list(map(lambda x : x.strip(), atomname))

    return atomname_list


def parse_xyz_for_scraps(xyz_fname) -> np.ndarray:
    """ parse coordinates  """

    x, y, z = np.loadtxt(xyz_fname, usecols=[1,2,3], skiprows=2, unpack=True)
    return x, y, z


def write_out_xyzfile(outfileName):
    """ To be extra sure that the correct xyz coordinates are parsed (in Aengstrom!), we do it ourselves.

        To not get in any trouble for whatever, the following script on the Orca-helpers gitlab is eerily close to what I wrote here.
        https://gitlab.gwdg.de/orca-helpers/orca-helpers/-/blob/master/orca_helpers/output_files/xyz_from_output/xyz_from_output.py

        Do with this information what you will, this function is all me.
        It is just string parsing from a file, so not all too surprising it is almost identical. """

    # Friendly neighbourhood message
    print("Parsing xyz coordinates from " + outfileName + "\n\n")


    # Parse the elements and coordinates from the MYFILE.out
    startParse = False
    xyzlist = list()
    with open(outfileName, "r") as OUTFILE :
        for line in OUTFILE:

            if line[:33] == "CARTESIAN COORDINATES (ANGSTROEM)" :
                next(OUTFILE)   # next() the iterator because it is a string not useful to us
                startParse = True
                continue

            if startParse:
                if len(line.split()) == 0 :
                    break
                else :
                    xyzlist.append(line)


    natoms = len(xyzlist)

    # write out elements and xyz coordinates to MYFILE_PARSED.xyz
    basename = outfileName.split(".")[0]
    with open("./" + basename + "_PARSED.xyz", "w") as xyzFILE:
        # amount of atoms
        xyzFILE.write(str(natoms) + "\n")
        xyzFILE.write("Coordinates from ORCA-job " + outfileName + "\n")

        for string in xyzlist:
            xyzFILE.write(string)

    print(f"Writing {basename}_PARSED.xyz \n\n")


def assign_vdW_radii(atomname_list : list) -> np.array :
    """ Radius data retrieved from antechamber and tleap.
        Converted prmtop and crd file, through pdb4amber, to a pqr format with atom radii.


        DATA : https://periodictable.com/Properties/A/VanDerWaalsRadius.v.html """

    radii_dict = {
            "C" : 1.700,
            "N" : 1.550,
            "O" : 1.520,
            "H" : 1.200,
            "P" : 1.800,
            }

    vdW_radii = np.zeros(shape=len(atomname_list))

    for i, atom in enumerate(atomname_list):
        vdW_radii[i] = radii_dict[atom[0]]

    return vdW_radii


def call_MSMS_for_grid(basename):
    """ Generate a couple of xyzr files to be read in by `msms`
        This will generated a couple *.vert files that will be used to parse the grid from. """
    # Parse atomname list from the pdb
    atomname_list = parse_pdb_for_scraps()
    # Parse radii by calling the pdb
    vdW_radii = assign_vdW_radii(atomname_list)
    # Parse xyz array
    x, y, z = parse_xyz_for_scraps(basename + "_PARSED.xyz")
    # Write an array from 1.4 to 2.0 to map electrostatic potential onto gridpoints, in A.U.
    factors = np.around(np.linspace(1.4, 2.0, num=4, endpoint=True), decimals=1)
    factors_name = [str(x).split(".")[0] + str(x).split(".")[1] for x in factors]

    # Start the shell script to run the `msms` program
    MSMS = open("./RUN_MSMS.sh", "w")
    MSMS.write( "#!/bin/bash\n"
                "executable=$(which msms)\n"
                "\n")

    # Write out the the files
    for i, factor in enumerate(factors):
        vdW_radii_sized = vdW_radii * factor
        with open("tmp_" + factors_name[i] + ".xyzr", "w") as tmpfile:
            for j in range(len(x)):
                tmpfile.write("{0:8.3f}{1:9.3f}{2:9.3f}{3:9.2f}\n".format(x[j], y[j], z[j], vdW_radii_sized[j]))

            # Write a line of command to the shell script
            MSMS.write(f"$executable -if tmp_{factors_name[i]}.xyzr -of tmp_{factors_name[i]} -probe_radius 1.4 -density 3.0 > tmp_{factors_name[i]}.out\n")
#            MSMS.write(f"$executable -if tmp_{factors_name[i]}.xyzr -of tmp_{factors_name[i]} -probe_radius 1.4 -density 1.0 > tmp_{factors_name[i]}.out\n")

    MSMS.close()
    print("Writing script to RUN_MSMS.sh\n")
    print("Writing xyzr file to tmp_*.vert\n")
    os.chmod(os.getcwd() + "/RUN_MSMS.sh", 0o775)


def prepare_vpot_input(basename) -> int:
    """ Parses all the *.vert files in the current working directory to make one giant grid.
        This is used as an input for orca_vpot (the xyz input). """

    # Needed to parse to tmp_*.vert files
    factors = np.around(np.linspace(1.4, 2.0, num=4, endpoint=True), decimals=1)
    tmp_fnames_prefix = [str(x).split(".")[0] + str(x).split(".")[1] for x in factors]

    sizeof = 0


    for i, prefix in enumerate(tmp_fnames_prefix):
        # Parse the array properly
        fname = "tmp_" + prefix + ".vert"
        tmp_array = np.loadtxt(fname, skiprows=3, dtype=float, usecols=[0,1,2])

        if i != 0:
            final_array = np.concatenate((final_array, tmp_array))
        elif i == 0:
            final_array = tmp_array



    ang_to_au = 1.0 / 0.5291772083
    x = np.array(final_array[:,0]) * ang_to_au
    y = np.array(final_array[:,1]) * ang_to_au
    z = np.array(final_array[:,2]) * ang_to_au

    sizeof = x.shape[0]

    with open("./" + basename + "_vpot.grid", "w") as GRID:
        GRID.write("{0:5d}\n".format(sizeof))

        for i, _ in enumerate(x):
            GRID.write("{0:12.6f}{1:12.6f}{2:12.6f}\n".format(x[i], y[i], z[i]))

    print("Writing "+ basename + "_vpot.grid\n")

    return sizeof


def write_ESP_grid_file(basename : str, npoints : int):
    """ Write the Electrostatic Potential Input file, required for `resp` """

    print("\n\nWriting esp.dat ... \n\n")


    # Parse xyz coordinates from the .xyz file, in Aengstrom
    x_atoms, y_atoms, z_atoms = parse_xyz_for_scraps("./" + basename + "_PARSED.xyz")

    # amount of atoms in the molecule
    natoms = len(x_atoms)

    ESP_INPUT = open("./esp.dat", "w")
    ESP_INPUT.write("{0:5d}{1:5d}{2:5d}\n".format(natoms, npoints, 0))

    ang_to_au = 1.0 / 0.5291772083
    # Convert to A.U.
    x_atoms *= ang_to_au
    y_atoms *= ang_to_au
    z_atoms *= ang_to_au

    # Read in orca_vpot's output (MYFILE_vpot.out)
    xgrid, ygrid, zgrid, vpot  = np.loadtxt("./" + basename + "_vpot.out", skiprows=1, unpack=True)

#    # Convert to A.U.
#    xgrid *= ang_to_au
#    ygrid *= ang_to_au
#    zgrid *= ang_to_au


    for i, _ in enumerate(x_atoms):
        ESP_INPUT.write("{0:32.7f}{1:16.7f}{2:16.7f}\n".format(x_atoms[i], y_atoms[i], z_atoms[i]))

    for j, _ in enumerate(vpot):
        ESP_INPUT.write("{0:16.7f}{1:16.7f}{2:16.7f}{3:16.7f}\n".format(vpot[j], xgrid[j], ygrid[j], zgrid[j]))

    ESP_INPUT.close()


def write_RESP_input_file(basename):
    """ Write the resp input file required to run the resp script """

    print("Writing to a standard resp.in file.\n"
            "   Check the file thoroughly before prompting it to resp.\n")

    # Initialise an atom dictionary to parse the correct amount of 
    atom_dict = {
            "H" : 1,
            "C" : 6,
            "N" : 7,
            "O" : 8,
            "F" : 9,
            "P" : 15,
            "S" : 16,
            "CL" : 17,
            "I" : 53,
            }

    # Parse correctly and remove any whitespace if there is
    # Collect the elements from the *_PARSED.xyz file
    elements = np.loadtxt(basename + "_PARSED.xyz", usecols=0, skiprows=2, dtype=str)
    element_list = list(map(lambda x : x.strip(), elements))

    natoms = len(element_list)

    # Start writing to the resp script
    RESP_INPUT = open("./resp.in", "w")
    RESP_INPUT.write(
        "TITLE\n"
        " &cntrl nmol=1, ihfree=1 \n"
        " /\n"
        "1.0\n"
        "TITLE - charge     natoms      \n"
            )

    RESP_INPUT.write("{0:5d}{1:5d}\n".format(0, natoms))

    for atom in element_list:
        if not atom in atom_dict:
            sys.exit(f"The following atom is not recognised by this script : {atom}. Either change the input or append it the existing dictionary ... (± line 262) ")
        RESP_INPUT.write("{0:5d}{1:5d}\n".format(atom_dict[atom], 0))

    RESP_INPUT.write("\n\n\n\n")
    RESP_INPUT.close()



def write_RESP_input_script():
    """ Write the shell script to prompt the command-line to resp
        This script assumes $AMBERHOME is a subdirectory of the $HOME path. Change if necessary.

        I could also parse the system's ~/.bashrc file for the $AMBERHOME, but I don't feel like doing that."""

    listHOME = os.listdir(os.path.expanduser("~"))

    AMBERHOME = [_dir for _dir in listHOME if "amber" in _dir]

    # If there are multiple instances of amber directories in home, just take the first one. Doesnt really matter which directory the resp script is called from
    if isinstance(AMBERHOME, list):
        AMBERHOME = os.path.expanduser("~") + "/" + AMBERHOME[0]

    RESP_SCRIPT = open("./resp_script.sh", "w")

    # Formalities
    RESP_SCRIPT.write(
        "#!/bin/bash\n"
        "\n"
        "source {}/amber.sh\n"
        "\n".format(AMBERHOME)
            )

    # commandline prompt
    RESP_SCRIPT.write(f"{AMBERHOME}/bin/resp -O \ \n"
                        "                    -i resp.in \ \n"
                        "                    -o resp.out \ \n"
                        "                    -p resp.dat \ \n"
                        "                    -e esp.dat \ \n"
                        "                    -t qout \ \n")
#                        "                    -q qin \n" )

    RESP_SCRIPT.close()

    # Modifiy permission to run as executable. Usually people do not like this to happen automatically, but I'm writing it anyway
    os.chmod(os.getcwd() + "/resp_script.sh", 0o775)

def REFERENCE_MICHEL_SANNER():
    """ MOST IMPORTANT PART IS THE FOLLOWING """

    reference = """
    AUTHOR OF THE MSMS SOFTWARE TO GENERATE MOLECULAR SURFACE VERTICES
    Michel F. Sanner,
    The Scripps Research Institute, La Jolla, California.


    MSMS is freely usable for academic and teaching purposes. For commercial research please contact Dr Sanner : sanner@scripps.edu

    If you use MSMS for your work please so kind and cite the following paper:
    Sanner, M. F., Olson A.J. & Spehner, J.-C. (1996). Reduced Surface: An Efficient Way to Compute Molecular Surfaces. Biopolymers 38:305-320.

    Link to the publication
    https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1097-0282(199603)38:3%3C305::AID-BIP4%3E3.0.CO;2-Y

    """

    print(reference)


def main():
    """ The steps are as follows :
        (1) Write out a MYFILE_PARSED.xyz file
        (2) Call the MSMS program to write out the vertices of the triangulated surface. This is your grid.
        (3) Call the orca_vpot program from orca to use the densities and the grid to calculate ESP charges for the gridpoints

        (4) Parse the data created and write out a esp.dat, resp.in and a shell script to run the resp script from AMBER

        (5) DO NOT FORGET TO QUOTE THE PAPER FROM PROF. DR. MICHEL SANNER """


    # Name of files in the current working directory
    try:
        basename = sys.argv[1]
    except IndexError:
        sys.exit("No file basename has been prompted as the first argument.\n"
                "Prompt the basename of the *.out file of the Hartree-Fock Single Point calculation.\n ")
    outfileName = basename + ".out"

    ### (1)
    # Initialise and check-up for writing out a proper MYFILE_PARSED.xyz in AENGSTROM
    # If the outfile is not in the current working directory, exit the script
    if os.path.isfile("./" + outfileName) :
        write_out_xyzfile(outfileName)
    else :
        sys.exit(f"Could not find the {basename}.out in the current working directory\n"
                "Make sure it is present in the current working directory for the script to work.\n")

    ### (2)
    # Parse everything up until actually running the RUN_MSMS.sh script.
    # NB : generates *_PARSED.xyz , `.vert` `.xyzr` files and a script to run the `msms` software
    call_MSMS_for_grid(basename)
    subprocess.check_call(["bash", "RUN_MSMS.sh"])
    print("Running RUN_MSMS.sh\n")

    # Prepare the orca_vpot input
    npoints = prepare_vpot_input(basename)

    ### (3)
    subprocess.check_call(['orca_vpot', basename + '.gbw', basename + '.scfp', basename + '_vpot.grid', basename + '_vpot.out'])

    ### (4)
    # Write everything to files formatted that are suitable for the `resp` software
    write_ESP_grid_file(basename, npoints)
    write_RESP_input_file(basename)
    write_RESP_input_script()

    ### (5)
    REFERENCE_MICHEL_SANNER()


#------------------------------------------------------------------------------------
#                                           MAIN
#------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

