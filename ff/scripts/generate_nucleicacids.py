import numpy as np
import sys, os, subprocess
import random


def check_dir(sequence):
    if not os.path.isdir(sequence):
        os.mkdir(sequence)
        return 1
    else:
        return


def write_to_NAB(sequence : str, fname) :

    filename = "./" + sequence + "/" + sequence + ".nab"
    with open(filename, "w") as gennuc:
        gennuc.write('molecule m;\n'
                    'm = fd_helix( "arna", "{}", "rna" );\n'
                    'putpdb( "{}/{}.pdb", m, "-wwpdb");\n'.format(sequence, sequence, sequence)
                    )

    fname.write("/home/jerome/amber18/bin/nab " + filename + "\n")
    fname.write("./a.out\n\n")


def get_random_string() :
    length = 12
    letters = ['a', 'g', 'c', 'u']

    result = "".join(random.choice(letters) for i in range(length))
    return result


def write_tleap(sequence : str, execute) :
    _tleap = open(sequence + '/tleap.in', 'w')
    _tleap.write('source leaprc.RNA.OL3\n'
                 'source leaprc.water.tip3p\n'
                 'complex = loadpdb ' + sequence + '.pdb\n'
                 'addions complex Na+ 0\n'
                 'addions complex Cl- 0\n'
                 'solvateoct complex TIP3PBOX 12\n'
                 'saveamberparm complex ' + sequence + '.prmtop ' + sequence + '.crd\n'
                 'savepdb complex leap.pdb\n'
                 'quit\n'
                 )
    _tleap.close()

    execute.write("$AMBERHOME/bin/tleap -f tleap.in\n\n")

def write_to_MD(sequence : str, execute) :

    ## MINIMISATION
    _min = open(sequence + '/min.in', 'w')
    _min.write('MINIMISATION RUN\n'
               '&cntrl\n'
               '  imin=1,\n'
               '  maxcyc=15000,\n'
               '  ncyc=7500,\n'
               '  ntb=1, ntc=2, ntf=2,\n'
               '  cut=10.0\n'
               '  ntpr=1000, ntwx=1000, ntwr=1000,\n'
               '    /\n'
               )
    _min.close()

    execute.write('$AMBERHOME/bin/pmemd.cuda -O -i min.in -o min.out -p ' + sequence + '.prmtop -c '
                    + sequence + '.crd -r min.ncrst -inf min.inf\n\n')

    ## HEATING
    _heat = open(sequence + '/heat.in', 'w')
    _heat.write('HEATING RUN\n'
                '&cntrl\n'
                '  imin=0, irest=0, ntx=1,\n'
                '  nstlim=75000, dt=0.002,\n'
                '  ntf=2, ntc=2,\n'
                '  cut=10.0, ntb=1,\n'
                '  vlimit=15.0,\n'
                '  ntt=3, gamma_ln=1.0,\n'
                '  tempi=0, temp0=300.0,\n'
                '  ntpr=1000, ntwx=1000, ntwr=1000,\n'
                '       /\n'
                )
    _heat.close()
    execute.write('$AMBERHOME/bin/pmemd.cuda -O -i heat.in -o heat.out -p ' + sequence + '.prmtop -c '
                  'min.ncrst -r heat.ncrst -x heat.nc -ref ' + sequence + '.crd -inf heat.inf\n\n')


    ## DENSITY
    _density = open(sequence + '/density.in', 'w')
    _density.write('DENSITY RUN\n'
                   '&cntrl\n'
                   '  imin=0, irest=1, ntx=5,\n'
                   '  nstlim=25000, dt=0.002,\n'
                   '  ntf=2, ntc=2,\n'
                   '  cut=10.0, ntb=2,\n'
                   '  ntp=1, taup=5.0,\n'
                   '  ntt=3, gamma_ln=2.0,\n'
                   '  tempi=300.0, temp0=300.0,\n'
                   '  ntpr=500, ntwx=500, ntwr=500,\n'
                   '       /\n'
                   )
    _density.close()
    execute.write('$AMBERHOME/bin/pmemd.cuda -O -i density.in -o density.out -p ' + sequence + '.prmtop -c '
                                      'heat.ncrst -r density.ncrst -x density.nc -ref ' + sequence + '.crd -inf density.inf\n\n')

    # EQUILIBRATION
    _equilibration = open(sequence + '/equi.in', 'w')
    _equilibration.write('EQUILIBRATION RUN\n'
                         '&cntrl\n'
                         '  imin=0, irest=1, ntx=5,\n'
                         '  nstlim=25000, dt=0.002,\n'
                         '  ntf=2, ntc=2,\n'
                         '  cut=10.0, ntb=2,\n'
                         '  ntp=1, taup=3.0,\n'
                         '  ntt=3, gamma_ln=2.0,\n'
                         '  tempi=300.0, temp0=300.0,\n'
                         '  ntpr=500, ntwx=500, ntwr=500,\n'
                         '       /\n'
                         )
    _equilibration.close()
    execute.write('$AMBERHOME/bin/pmemd.cuda -O -i equi.in -o equi.out -p ' + sequence + '.prmtop -c '
                                      'density.ncrst -r equi.ncrst -x equi.nc -ref ' + sequence + '.crd -inf equi.inf\n\n')

    # PRODUCTION
    prod_runs = 40              # DETERMINES THE AMOUNT OF NANOSECONDS THAT THE PRODUCTION RUN IS CARRIED OUT

    _production = open(sequence + '/prod.in', 'w')
    _production.write('PRODUCTION RUN\n'
                      '&cntrl\n'
                      '  imin=0, irest=1, ntx=5,\n'
                      '  nstlim=500000, dt=0.002, pencut=-0.001,\n'
                      '  ntf=2, ntc=2,\n'
                      '  cut=10.0, ntb=2,\n'
                      '  ntp=1, taup=2.0,\n'
                      '  ntt=3, gamma_ln=2.0,\n'
                      '  tempi=300.0, temp0=300.0,\n'
                      '  ntpr=5000, ntwx=5000, ntwr=20000,\n'
                      '       /\n'
                      )
    _production.close()

    executefile.write('$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod1.out -p ' + sequence + '.prmtop -c equi.ncrst -r prod1.ncrst -x prod1.nc '
                      '-ref ' + sequence + '.crd -inf prod.inf\n')

    for i in range(2, prod_runs + 1):
        j = i - 1
        execute.write('$AMBERHOME/bin/pmemd.cuda -O -i prod.in -o prod' + str(i) + '.out -p ' + sequence + '.prmtop '
                        '-c prod' + str(j) + '.ncrst -r prod' + str(i) + '.ncrst -x prod' + str(i) + '.nc '
                        '-ref ' + sequence + '_complex.crd -inf prod.inf\n')


if __name__ == "__main__":
    # initialise runfile

    genDuplex = open("./generateDuplexes.sh", "w")
    # LOOP
    # The range generates X amount of double strands
    for i in range(10):

        # Generate duplex
        duplex = get_random_string()

        # Check if the directory already exists
        if not check_dir(duplex):
            continue

        # Create sourcefile, in bash, to be able to source Amber
        sourcefile = open("./" + duplex + "/sourcfile.sh", "w")
        sourcefile.write("#!/bin/bash\n")
        sourcefile.write('export AMBERHOME="$HOME"/amber22\n')
        sourcefile.write("source $AMBERHOME/amber.sh\n\n\n")

        # write nab file to the directory, that generates the duplex
        write_to_NAB(duplex, sourcefile)

        # Create runfile to run the MD commands for the simulation
        executefile = open("./" + duplex + "/run_MD.sh", "w")
        executefile.write("#!/bin/bash\n")
        sourcefile.write('export AMBERHOME="$HOME"/amber22\n')
        executefile.write("source $AMBERHOME/amber.sh\n\n\n")

        write_tleap(duplex, executefile)
        write_to_MD(duplex, executefile)

        genDuplex.write("echo '" + duplex + "' | tee -a progress.log\n")
        genDuplex.write("(time  $AMBERHOME/bin/nab " + "./" + duplex + "/" + duplex + ".nab) 2>&1 | tee -a progress.log\n")
        genDuplex.write("(time  ./a.out) 2>&1 | tee -a progress.log \n")
        genDuplex.write("echo '---------------' >> progress.log\n")

        sourcefile.close() ; executefile.close

    genDuplex.close()
    subprocess.run(["bash", "./generateDuplexes.sh"])
