
from fileinput import input as input_
from glob import glob as glob_
import numpy as np
import pandas as pd


# Globular loadtxt
def loadtxt(fname, dtype=float, glob=False, verbose=False):
    """
    A refactoring of :ref:`numpy.loadtxt` that allows for globbing files.

    Parameters
    ----------
    fname : file, str, or pathlib.Path
        Name of file.
    dtype : str or object
        File type.
    glob : bool
        Does `fname` need to be globbed?
    verbose : bool
        Should information about the read-in be displayed?

    Returns
    -------
    pandas.Series
        Read file
    """

    # If glob, change fname to include all globbed files
    if glob:
        # Glob first; if glob is empty, throw an error
        fname_glob = glob_(fname)
        if not fname_glob:
            raise FileNotFoundError(fname)

        # Sort glob
        fname_glob = sorted(fname_glob)

        # Output if verbose
        if verbose:
            print(f'glob: {list(fname_glob)}')

        # Update fname to include all globbed files
        fname = input_(fname_glob)

    # Utilize numpy to read-in the file(s)
    data = np.loadtxt(fname, dtype=dtype)

    # If verbose, note the shape of the data
    if verbose:
        print(f'file loaded with shape {data.shape}')

    # Return
    return data


def read_table(fname, glob=False, verbose=False, **kwargs):
    """
    Read table into :class:`pandas.DataFrame`.

    Parameters
    ----------
    fname
    glob
    index_col
    verbose

    Returns
    -------

    """

    # If glob, change fname to include all globbed files
    if glob:
        # Glob first; if glob is empty, throw an error
        fname_glob = glob_(fname)
        if not fname_glob:
            raise FileNotFoundError(fname)

        # Sort glob
        fnames = sorted(fname_glob)

        # Output if verbose
        if verbose:
            print(f'glob: {fnames}')

    # Otherwise, turn fname into a list
    # TODO evaluate if creating this list is right, or if we should short-circuit the read-in
    else:
        fnames = [fname]

    # Cycle over fnames and read in
    dataframes = [pd.read_table(_, **kwargs) for _ in fnames]

    # Concatenate
    dataframe = dataframes[0] if len(dataframes) == 1 else pd.concat(dataframes)

    # Return
    return dataframe
