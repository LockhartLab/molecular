"""
transform.py

language: Python
version: 3.x
author: C. Lockhart <chris@lockhartlab.org>
"""

import logging
from numba import njit
import numpy as np

# Get the molecular.transform logger
logger = logging.getLogger('molecular.transform')


# Compute the center of the Trajectory
def center(a, weights=None):
    r"""
    Compute the center `m` of Trajectory `a`.

    .. math :: m = \frac{\sum_{i=1}^N w_i (x_i + y_i + z_i)} {\sum_{i=1}^N w_i}

    Here, :math:`w_i` refers to weights applied to individual atoms :moath:`i`.

    Parameters
    ----------
    a : molecular.Trajectory
    weights : array-like
        (Optional) Weights for atoms in the trajectory. Follows the definition from :ref:`numpy.average`.

    Returns
    -------
    numpy.ndarray
        Center of the Trajectory.
    """

    # Compute center
    m = np.average(a.xyz, weights=weights, axis=1, returned=False)

    # Update log
    logging.info(f'computed center of {a.designator}')

    # Return
    return m


def fit(a, b):
    """
    Fit Trajectory `b` to Trajectory `a` in a pairwise fashion.

    Parameters
    ----------
    a : molecular.Trajectory
    b : molecular.Trajectory

    Returns
    -------

    """

    # Number of structures
    n_structures = a.n_structures
    if n_structures != b.n_structures:
        raise AttributeError('a and b must have the same number of structures')

    # Number of atoms
    n_atoms = a.n_atoms
    if n_atoms != b.n_atoms:
        raise AttributeError('a and b must have the same number of atoms')

    # Number of dimensions
    n_dim = 3

    # Move structures to center
    a_xyz_center = center(a)
    b_xyz_center = center(b)
    a_xyz = a.xyz - a_xyz_center
    b_xyz = b.xyz - b_xyz_center

    # Transpose b
    b_xyz_transpose = np.transpose(b_xyz, axes=[0, 2, 1])

    # Compute covariance matrix
    covariance_matrix = np.matmul(b_xyz_transpose, a_xyz)

    # Get rotation matrix
    rotation_matrix = _get_rotation_matrix(covariance_matrix)

    # Perform optimal rotation
    b_xyz = np.matmul(b_xyz, rotation_matrix)

    # Move `b` to `a` center
    b_xyz = b_xyz + a.center()

    # Return
    return b_xyz


# Center Trajectory
# def center(a, weights=None, inplace=False):
#     r"""
#     Compute the center of Trajectory `a`.
#
#     .. math :: center = \frac{1}{N} \Epsilon w_i (x_i + y_i + z_i)
#
#     Parameters
#     ----------
#     a : molecular.Trajectory
#     weights : array-like
#         (Optional) Weights for atoms in the trajectory. Follows the definition from :ref:`numpy.average`.
#     inplace : bool
#         Should the Trajectory be centered in-place? (Default: False)
#
#     Returns
#     -------
#     molecular.Trajectory
#         Trajectory centered at origin.
#     """
#
#     # Compute center
#     xyz = np.average(a.xyz, weights=weights, axis=1, returned=False)
#
#     # Update log
#     logging.info(f'computed center of {a.designator}')
#
#     # Copy if not an inplace transformation
#     if not inplace:
#         a = a.copy()
#
#     # Update the coordinates
#     a.xyz = xyz
#
#     # Return the copy if necessary
#     if not inplace:
#         return a
#


# Move Trajectory
def move(a, by=None, to=None, return_copy=False):
    r"""
    Move Trajectory `a` either `by` a constant or `to` a location. The two arguments `by` and `to` cannot be set
    simultaneously.

    Parameters
    ----------
    a : molecular.Trajectory
    by : tuple
        (Optional) Move `a` by a constant.
    to : tuple
        (Optional) Move `a` to a location.
    return_copy : bool
        Should we perform the translation on a copy of the Trajectory (Default: False)

    Returns
    -------
    None or molecular.Trajectory
        If `return_copy` is set, we return the moved Trajectory.
    """

    # If both `by` and `to` are None (or both set), we don't know what to do
    if (by is None and to is None) or (by is not None and to is not None):
        raise AttributeError('either `by` or `to` must be set')

    # If `to` is set, let's refactor this as `by`
    if to is not None:
        m = center(a)
        by = to - m  # xyz - center + to is the idea

    # `by` must be the same shape as n_dim
    # if by.ndim != 1 or by.shape[0] != a.n_dim:
    #     raise AttributeError(f'by must contain {a.n_dim} elements')

    # Should we make a copy?
    if return_copy:
        a = a.copy()

    # Perform translation
    # TODO need to come up with robust tests about when array links are preserved
    # noinspection PyProtectedMember
    a._xyz[:] = a._xyz[:] + by  # Preserves the transformation to any linked objects

    # Update log
    logging.info(f'translated {a.designator} by {by}')

    # Return if necessary
    if return_copy:
        return a



def _fit(a_xyz, b_xyz, backend='python'):
    """
    `a_xyz` and `b_xyz` must be paired.

    Parameters
    ----------
    a_xyz : numpy.ndarray
    b_xyz : numpy.ndarray
    backend : str

    Returns
    -------

    """

    # Extract dimensions
    n_structures, n_atoms, n_dim = a_xyz.shape

    # Sanity check
    if n_structures != b_xyz.shape[0]:
        raise AttributeError('a_xyz and b_xyz must have the same number of structures')
    if n_atoms != b_xyz.shape[1]:
        raise AttributeError('a_xyz and b_xyz must have the same number of atoms')

    # Centers
    a_xyz_center = np.tile(a_xyz.mean(axis=1), n_atoms).reshape(n_structures, n_atoms, n_dim)
    b_xyz_center = np.tile(b_xyz.mean(axis=1), n_atoms).reshape(n_structures, n_atoms, n_dim)

    # Move structures to center
    a_xyz = a_xyz - a_xyz_center
    b_xyz = b_xyz - b_xyz_center

    # Transpose b
    b_xyz_transpose = np.transpose(b_xyz, axes=[0, 2, 1])

    # Compute covariance matrix
    covariance_matrix = np.matmul(b_xyz_transpose, a_xyz)

    # Get rotation matrix
    rotation_matrix = _get_rotation_matrix(covariance_matrix)

    # Perform optimal rotation
    b_xyz = np.matmul(b_xyz, rotation_matrix)

    # Move b to a center
    b_xyz = b_xyz + a_xyz_center

    # Return
    return b_xyz


@njit
def _get_rotation_matrix(covariance_matrix):
    """
    Compute the rotation matrix from the method from Kabsch (1976) Acta Cryst. A32.



    Parameters
    ----------
    covariance_matrix

    Returns
    -------

    """
    n_structures = covariance_matrix.shape[0]
    rotation_matrix = []
    for i in range(n_structures):
        u, s, v = np.linalg.svd(covariance_matrix[i, :, :], full_matrices=True)
        s = np.diag(np.array([1., 1., np.sign(np.linalg.det(np.dot(u, v)))]))
        r = np.dot(np.dot(u, s), v)
        rotation_matrix.append(r)
    return rotation_matrix


if __name__ == '__main__':
    import molecular as mol

    trj = mol.read_pdb('../tests/samples/trajectory.pdb')
    trj.to_center(inplace=True)
    assert (trj.center() < 1e-7).max()
