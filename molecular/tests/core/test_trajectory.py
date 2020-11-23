"""
test_trajectory.py

language: Python3
author: C. Lockhart <chris@lockhartlab.org>
"""

# from molecular.core import Structure, Trajectory
from molecular.core import Trajectory
from molecular.io import read_pdb

import numpy as np
import os.path
from pytest import fixture


# Path to samples
samples = os.path.join('molecular', 'tests', 'samples')


@fixture
def load_pdb():
    return False


def test_recenter():
    pass
