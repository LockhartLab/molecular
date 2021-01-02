
from molecular.analysis import order
from molecular.misc import experimental


# TODO compute with Voronoi tesselation
# how to remove protein, other components?
@experimental
def area_per_lipid():
    pass


# TODO compute carbon-deuterium order parameter -- perhaps generalize this?
@experimental
def scd(carbon, hydrogen):
    r"""

    What form should this be in?  # str by # lipids by lipid length?? Some lipids will have unequal length...

    [ [hydrogen1 - carbon, hydrogen2 - carbon] ]



    ..math :: S_{CD} = \frac{3}{2}cos^2 \theta - \frac{1}{2}

    Returns
    -------

    """

    pass


@experimental
def tilt():
    pass
