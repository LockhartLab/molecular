
import numpy as np
from scipy.io import FortranFile

cimport numpy as np

def _read_dcd(fname):
    # Open FortranFile buffer
    buffer = FortranFile(fname, 'r')

    # Header
    cdef int n_str, fixed
    header, n_str, _, fixed, _ = buffer.read_record('4a', 'i', '7i', 'i', '11i')
    if header[0] != b'CORD' or fixed[0] != 0:
        raise IOError('cannot parse DCD file')
    n_str = n_str[0]

    # Title
    cdef int n_titles
    n_titles, titles = buffer.read_record('i', '2a80')
    if n_titles != 2:
        raise IOError('cannot parse DCD file')

    # Atoms
    cdef n_atoms
    n_atoms = buffer.read_record('i')
    n_atoms = n_atoms[0]

    # Box and coordinate information
    cdef np.ndarray box = np.zeros((n_str, 3))
    cdef np.ndarray xyz = np.zeros((n_str, n_atoms, 3))
    for i in range(n_str):
        box[i, 0], _, box[i, 1], _, box[i, 2] = buffer.read_record('f8', 'f8', 'f8', '2f8', 'f8')
        xyz[i, :, 0] = buffer.read_record(f'{n_atoms}f')
        xyz[i, :, 1] = buffer.read_record(f'{n_atoms}f')
        xyz[i, :, 2] = buffer.read_record(f'{n_atoms}f')

    # Close out buffer
    buffer.close()

    # Return
    return box, xyz

