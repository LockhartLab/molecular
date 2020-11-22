"""
test_energy.py
written in Python3
author: C. Lockhart <chris@lockhartlab.org>
"""


from molecular.energy import *

from hypothesis import given
from hypothesis import strategies as st
import numpy as np
from unittest import TestCase


class TestEnergy(TestCase):
    @given(st.integers(min_value=1, max_value=1000),
           st.integers(min_value=1, max_value=3),
           st.floats(min_value=0., max_value=100.),
           st.floats(min_value=0., max_value=100.))
    def test_bond(self, n_atoms, n_dim, r0, k):
        # Create points
        a = np.random.rand(n_atoms, n_dim)
        b = np.random.rand(n_atoms, n_dim)
        c = np.random.rand(n_atoms, n_dim)
        d = np.random.rand(n_atoms, n_dim)

        # Get vectors
        u = b - a
        v = c - b
        w = d - c

        # Compute unit vectors
        u_norm = u / _norm(u, n_atoms)
        v_norm = v / _norm(v, n_atoms)

        # Sanity
        r = np.sqrt(np.sum(np.square(b - a), axis=1))
        x = (b - a) / r[:, None]
        y = np.array([(b - a)[:, i] / r for i in range(a.shape[1])]).T
        np.testing.assert_array_equal(x.ravel(), np.array(y).ravel())

        # Compute bond metrics
        vector = b - a
        distance = np.sqrt(np.sum(np.square(vector), axis=1))
        offset = distance - r0
        energy = 0.5 * k * offset ** 2
        force = -k * offset[:, None] * vector / distance[:, None]
        term = Bond(a, b, r0, k)

        # Test bond
        np.testing.assert_array_almost_equal(vector, term._vector)
        np.testing.assert_array_almost_equal(distance, term._distance)
        np.testing.assert_array_almost_equal(offset, term._offset)
        np.testing.assert_array_almost_equal(energy, term.energy)
        np.testing.assert_array_almost_equal(force, term.force)

        # Angle
        if n_dim >= 2:
            angle = np.arccos(_dot(u_norm, v_norm, n_atoms))
            offset = angle - r0
            energy = 0.5 * k * np.square(offset)
            # force = -k * offset[:, None] * vector / distance[:, None]
            term = Angle(a, b, c, r0, k)
            np.testing.assert_array_almost_equal(angle, term._angle)
            np.testing.assert_array_almost_equal(offset, term._offset)
            np.testing.assert_array_almost_equal(energy, term.energy)


def _dot(u, v, n_points):
    if n_points == 0:
        result = np.vdot(u, v)
    else:
        result = np.sum(np.multiply(u, v), axis=1)
    return result


def _norm(u, n_points):
    if n_points == 0:
        result = np.linalg.norm(u)
    elif n_points == 1:
        result = np.linalg.norm(u, axis=1)
    else:
        result = np.linalg.norm(u, axis=1).reshape(-1, 1)
    return result
