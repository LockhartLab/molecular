
import sys
sys.path.insert(0, '..')
from molecular import read_pdb

pdb = read_pdb('../samples/trajectory.pdb')
print(pdb.n_structures)
