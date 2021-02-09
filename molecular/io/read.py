
"""
read.py
written in Python3
author: C. Lockhart <chris@lockhartlab.org>
"""

from molecular.core import Topology, Trajectory

# from molecular.io.fortran.read_dcd import read_dcd as _read_dcd  # noqa

import numpy as np
# from numpy.lib.recfunctions import drop_fields, structured_to_unstructured
import pandas as pd
import re
from scipy.io import FortranFile
from typelike import ArrayLike


# Read PDB
# TODO currently the only backend will by pandas; in future, expand to Cython or C or Fortran backend
def read_pdb(fname, backend='python'):
    """
    Read PDB file and return Trajectory

    PDB format can be described in depth at `<http://www.wwpdb.org/documentation/file-format>`_.

    Parameters
    ----------
    fname : str
        Name of PDB file to be read
    backend : str
        (Default: 'pandas')

    Returns
    -------
    molecular.Trajectory
        Trajectory of PDB
    """

    # Make sure we know we're using the pandas backend
    if backend.lower() != 'python':
        raise AttributeError('only python backend presently supported')

    # Open file, read in all records
    with open(fname, 'r') as stream:
        records = stream.read()
        # records = _atom_reader(buffer)

    # Return
    return _read_pdb(records)


def read_peptide_sequence(residues):
    """
    Create trajectory from an amino acid sequence.

    Parameters
    ----------
    residues : str or ArrayLike
        List of amino acids.

    Returns
    -------
    molecular.Trajectory
        Instance of Trajectory.
    """

    # Convert to list if not already list
    if not isinstance(residues, ArrayLike):
        residues = list(residues)

    # Make sure we're in the right format
    letter_to_code = {
        'A': 'ALA',
        'R': 'ARG',
        'N': 'ASN',
        'D': 'ASP',
        'C': 'CYS',
        'Q': 'GLN',
        'E': 'GLU',
        'G': 'GLY',
        'H': 'HIS',
        'I': 'ILE',
        'L': 'LEU',
        'K': 'LYS',
        'M': 'MET',
        'F': 'PHE',
        'P': 'PRO',
        'S': 'SER',
        'T': 'THR',
        'W': 'TRP',
        'Y': 'TYR',
        'V': 'VAL'
    }
    residues = [letter_to_code.get(residue, residue) for residue in residues]

    # Create Trajectory
    topology = Topology()
    topology._data = pd.DataFrame({'residue_id': np.arange(len(residues)), 'residue': residues})

    # Create Trajectory
    trajectory = Trajectory(xyz=np.zeros((1, len(residues), 3)), topology=topology)

    # Return
    return trajectory


def _read_pdb(records):
    # Filter out CRYST1 records
    cryst = re.sub(r'^(?!CRYST1).*$', '', records, flags=re.MULTILINE).lstrip()
    n_structures = len(re.split('\WCRYST1', cryst))

    # Filter out atom records
    # TODO this will be slow for large PDB files; perhaps move to Cython or C backend
    atoms = re.sub(r'^(?!ATOM).*$', '', records, flags=re.MULTILINE).replace('\n\n', '\n').lstrip()

    # Sections of PDB
    sections = np.array([
        (6, 'record', '<U6'),
        (5, 'atom_id', 'int'),
        (5, 'atom', '<U5'),
        (5, 'residue', '<U5'),
        (1, 'chain', '<U1'),
        (4, 'residue_id', 'int'),
        (4, 'blank', '<U1'),
        (8, 'x', 'float'),
        (8, 'y', 'float'),
        (8, 'z', 'float'),
        (6, 'alpha', 'float'),
        (6, 'beta', 'float'),
        (10, 'segment', '<U9'),
        (2, 'element', '<U2')
    ], dtype=[('width', 'i1'), ('column', '<U10'), ('type', '<U10')])

    # Read in
    data = np.genfromtxt(atoms.split('\n'), delimiter=sections['width'], dtype=sections['type'], autostrip=True)
    # data.dtype.names = sections['column']
    data = pd.DataFrame(data.tolist(), columns=sections['column'])

    # Drop extraneous columns
    # data = drop_fields(data, 'blank')
    data = data.drop(columns='blank')

    # TODO this should also be done for residue_id probably
    # If the PDB starts at atom_id = 1, change to 0-index
    # if data['atom_id'].min() == 1:
    #     data['atom_id'] -= 1

    # Determine number of structures in PDB
    # atom_id = data['atom_id'].values
    # pos_atom_id = atom_id[atom_id >= 0]
    # _, atom_counts = np.unique(pos_atom_id, return_counts=True)
    # n_structures = np.unique(atom_counts)
    # # n_structures = data.pivot_table(index='atom_id', values='record', aggfunc='count')['record'].unique()
    # if len(n_structures) != 1:
    #     raise AttributeError('inconsistent record counts in PDB, %s' % n_structures)
    #     # raise AttributeError('inconsistent record counts in PDB')
    # n_structures = n_structures[0]

    # Separate out dynamic columns for Trajectory and static Topology data
    dynamical_columns = ['x', 'y', 'z']
    # static_columns = [column for column in data.dtype.names if column not in dynamical_columns]
    static_columns = [column for column in data.columns if column not in dynamical_columns]

    # Renumber atom_id
    if np.mod(len(data), n_structures) != 0:
        raise AttributeError('len(data) must be evenly divisible by n_structures, (%s, %s)' % (len(data), n_structures))
    # TODO, why does this index start at 1? I changed this 1/28/21 but might break something
    # data['atom_id'] = np.array(np.tile(np.arange(1, len(data) / n_structures + 1), n_structures), dtype='int')
    data['atom_id'] = np.array(np.tile(np.arange(len(data) / n_structures), n_structures), dtype='int')

    # Create Topology first
    # TODO what happens when alpha and beta differ by structures? Should these be stored in Trajectory?
    data['alpha'] = 0.
    data['beta'] = 0.
    # topology = Topology(np.unique(data[static_columns]))
    topology = Topology(data[static_columns].drop_duplicates())

    # Next create Trajectory (the result)
    # n_atoms = len(np.unique(data['atom_id'].values))
    n_atoms = data['atom_id'].nunique()
    # result = Trajectory(structured_to_unstructured(data[dynamical_columns]).reshape(n_structures, n_atoms, 3),
    #                     topology=topology)
    result = Trajectory(data[dynamical_columns].values.reshape(n_structures, n_atoms, 3), topology=topology)

    # Return
    return result


# @jit(nopython=False)
# def _atom_reader(buffer):
#     result = ''
#     for line in buffer.readlines():
#         if line[:4] == 'ATOM':
#             result += line
#     return result

# Read DCD
# could https://www.pytables.org/usersguide/tutorials.html speed this up?
def read_dcd(fname, topology=None, backend='scipy'):
    """
    Read in DCD file with `fname`. This function is partially based off
    https://www.ks.uiuc.edu/Research/namd/wiki/index.cgi?ReadingDCDinFortran
    
    Parameters
    ----------
    fname : str
        Name of DCD file.
    topology : Topology
        (Optional) Topology file to load for coordinates.
    backend : str
        How to load the DCD? (Default: "scipy")

    Returns
    -------
    Trajectory
        Instance of Trajectory object.
    """

    # Convert backend to lowercase
    backend = backend.lower()

    # Our home-grown Scipy backend
    if backend in 'scipy':
        # Open FortranFile buffer
        buffer = FortranFile(fname, 'r')

        # Header
        header, n_str, _, fixed, _ = buffer.read_record('4a', 'i', '7i', 'i', '11i')
        if header[0] != b'CORD' or fixed[0] != 0:
            raise IOError('cannot parse DCD file')
        n_str = n_str[0]

        # Title
        n_titles, titles = buffer.read_record('i', '2a80')
        if n_titles != 2:
            raise IOError('cannot parse DCD file')

        # Atoms
        n_atoms = buffer.read_record('i')
        n_atoms = n_atoms[0]

        # Box and coordinate information
        box = np.zeros((n_str, 3))
        xyz = np.zeros((n_str, n_atoms, 3))
        for i in range(n_str):
            box[i, 0], _, box[i, 1], _, box[i, 2] = buffer.read_record('f8', 'f8', 'f8', '2f8', 'f8')
            xyz[i, :, 0] = buffer.read_record(f'{n_atoms}f')
            xyz[i, :, 1] = buffer.read_record(f'{n_atoms}f')
            xyz[i, :, 2] = buffer.read_record(f'{n_atoms}f')

        # Close out buffer
        buffer.close()

    # MDAnalysis for comparison and error checking
    elif backend in 'mdanalysis':
        from MDAnalysis import Universe  # noqa

        # MDAnalysis has a Cython backend
        universe = Universe(fname, in_memory=True)
        box = universe.trajectory.dimensions_array[:, :3]
        xyz = universe.trajectory.coordinate_array

    # Unknown
    else:
        raise AttributeError(f'unknown backend {backend}')

    # Build Trajectory
    return Trajectory(xyz, box=box, topology=topology)
