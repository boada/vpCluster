from numpy.random import choice
from numpy import sort

def bootstrap(data, statistic, resamples=1000, alpha=0.05, **args):
    """ Returns the bootstrap estimate of the confidence interval for the given
    statistic. The confidence interval is given by 100*(1-alpha). Passes a 1d
    array to the function, statistic. Any arguments needed by statistic are
    passed by **args.

    @type data: list
    @param data: The data on which the given statistic is calculated
    @type statistic: function
    @param statistic: The statistic desired
    @type resamples: int
    @param resamples: The number of bootstrap resamplings
    @type alpha: float
    @param alpha: The confidence interval given by 100*(1-alpha), 95% default
    @type **args: Keywords
    @param **args: Arguments needed by the statistic function
    @rtype: tuple
    @return: (Lower Interval, Upper Interval)

    """

    samples = choice(data, size=(resamples, len(data)), replace=True)
    stat = sort([statistic(row, **args) for row in samples])
    return (stat[int((alpha/2.0) * resamples)],
            stat[int((1-alpha/2.0) * resamples)])

