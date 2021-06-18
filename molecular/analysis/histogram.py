
import numpy as np


class Histogram:
    def __init__(self, bins, counts, weights=None):
        self._bins = bins
        self._counts = counts
        # FIXME for non-ihist, the width will impact the mass function
        if weights is None:
            weights = 1.
        self._probs = counts / np.sum(counts)  # pmf

    def plot(self):
        import uplot as u
        u.plot(self.to_frame())

    # FIXME when count = False and probability = False, this will break because there will be a pandas index with no
    #  other content.
    def to_frame(self, count=True, probability=False):
        import pandas as pd

        df = pd.DataFrame({'bin': self._bins})
        if count:
            df['count'] = self._counts
        if probability:
            df['probability'] = self._probs

        return df.set_index('bin')


def hist(a, bins):
    """
    Classical histogram, numpy style.

    Parameters
    ----------
    a
    bins

    Returns
    -------

    """

    pass


def ihist(a, width, right=False):
    """
    Create an interval histogram, i.e., given a known `width`, count all the occurrences of `a` in all bins :math:`[
    width, 2*width)`. By default, `right = False`. If true, then we look at :math:`(width, 2*width]` instead.

    Parameters
    ----------
    a : 1D array
        Collection of values to turn into histogram.
    width : number
        The interval that should be binned.
    right : bool
        Left or right inclusive? If False, left inclusive. (Default: False)

    Returns
    -------
    Histogram
        Object that contains the histogram.
    """

    # Take `a` and get binned values
    func = np.floor if not right else np.ceil
    a = np.copy(func(a / width) * width)

    # Create Histogram object and return
    return Histogram(*np.unique(a, return_counts=True))


if __name__ == '__main__':
    a = np.random.rand(100000)
    print(ihist(a, 0.1).to_pandas(probability=True, count=False))