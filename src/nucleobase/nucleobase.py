# first file of the WIP
from nucleobase.utils_nucleobase import Nucleobase



def modify_nucleobases(PDB_NBASE_FNAME: str, LIST_OF_NBASE_MODIFICATIONS : list[Nucleobase]) -> None : 
    """ 
    take in the CLI inputs
    -> pdb file 
        
    -> pre-processed CLI inputs as a Class
        --nucleobase 
        -position 15.A  (required, where chain is optional)
        -mod Psi        (required)
        -resname pU     (optional, takes on resname of the current residue in the original pdb)
        -orientation    (optional)


    return : Outputs a new pdb file 
    """



