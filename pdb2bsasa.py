"""
Given a .pdb file as input, compute the Buried Solvent Accessible Surface Area (BSASA)
for all chain pairs, and write to stdout. 

BSASA is the sum of monomer SASA to binary complex SASA.
BSASA = 0.5 * (SASA_A + SASA_B - SASA_AB)

SASA is calculated with default parameters.
- 1.4 Angstrom probe radius
- Protor Radii set 

"""
import biotite.structure.io as strucio
import biotite.structure as struc
from itertools import combinations
import math
import numpy as np
import sys
import click

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs) 

@click.command()
@click.option("--pdb-path", type=str)
@click.option("--pairwise-only", is_flag=True, default=False, help="enable for binary complexes only")
@click.option("--sep", type=str, default=",", help="seperator character")
@click.option("--header", is_flag=True, default=False)
def main(pdb_path, pairwise_only, sep, header):
    if header:
        header = f"Path{sep}ChainA{sep}ChainB{sep}BSASA{sep}SASA_A{sep}SASA_B{sep}SASA_AB"
        print(header)
    if pdb_path:
        atom_array = strucio.load_structure(pdb_path)
        chain_set = set(atom_array.chain_id)
        nchains = len(chain_set)
        if pairwise_only:
            assert nchains == 2, f"Expected a binary complex, found {nchains} in {pdb_path}"
        for (A, B) in combinations(chain_set, 2):
    
            chain_A = atom_array[atom_array.chain_id == A]
            chain_B = atom_array[atom_array.chain_id == B]
            #Per atom SASA
            chain_A_atom_sasa = struc.sasa(chain_A) # Protor default radii
            chain_B_atom_sasa = struc.sasa(chain_B)
    
            chain_A_sasa = struc.apply_chain_wise(chain_A, chain_A_atom_sasa, np.sum)
            chain_B_sasa = struc.apply_chain_wise(chain_B, chain_B_atom_sasa, np.sum)
            sel = (atom_array.chain_id == A) | (atom_array.chain_id == B) 
            chain_AB = atom_array[sel] 
            chain_AB_atom_sasa = struc.sasa(chain_AB)
            chain_AB_sasa = np.sum(chain_AB_atom_sasa)
            bsasa = (chain_A_sasa + chain_B_sasa - chain_AB_sasa)
            out = f"{pdb_path}{sep}{A}{sep}{B}{sep}{bsasa.item()}{sep}{chain_A_sasa.item()}{sep}{chain_B_sasa.item()}{sep}{chain_AB_sasa.item()}"
            print(out)
        
if __name__ == "__main__":
    main()
