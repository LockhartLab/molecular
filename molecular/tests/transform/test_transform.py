import molecular as mol
import os

# Path to samples
if __name__ == '__main__':
    samples = os.path.abspath(os.path.join('..', 'samples'))
else:
    samples = os.path.abspath(os.path.join(__file__, '..', '..', 'samples'))


def load_pdb():
    return mol.read_pdb(os.path.join(samples, 'trajectory.pdb'))


def test_center():
    trj = load_pdb()

    trj.move(to=(0, 0, 0))
    assert (trj.center() < 1e-7).max()
    assert (mol.center(trj) < 1e-7).max()


def test_move():
    trj = load_pdb()

    trj.move(to=(0, 0, 0))
    assert (trj.center() < 1e-7).max()

    mol.move(trj, to=(0, 0, 0))
    assert (trj.center() < 1e-7).max()

    trj.move(by=(5, 5, 5))
    assert (trj.center() == 5).max()

    mol.move(trj, by=(0, 0, 0))
    assert (trj.center() == 5).max()
