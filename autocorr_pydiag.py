#!/usr/bin/env python3
import numpy as np

## stolen from gene3d-dev/diag-python/pydiag/utils/averages.py 
def mytrapz(yvar, timefld):
    """ Trapezoid rule, adjusted for single timestep

    Operates on the first dimension of yvar (typically time)
    hence the name timefld for the integration variable samples
    :param yvar: Variable to integrate over
    :param timefld: The integration variable (time)
    :returns: The integration result (dim(yvar)-1)
    """
    timefld = np.array(timefld)
    yvar = np.real_if_close(yvar)
    if timefld.size == 1:
        return np.atleast_2d(yvar)[0]
    else:
        if yvar.shape[0] != len(timefld):
            raise ValueError("First dimension of yvar and timefld do not match")
        tmpt = timefld/(timefld[-1] - timefld[0])
        return np.trapz(yvar, x=tmpt, axis=0)


## stolen from gene3d-dev/diag-python/pydiag/utils/errors.py
def autocorrtime_1d(data, timefld):
    """ Correlation time calculation for 1d arrays (only time)

    :param data: 1d array of the data to correlate
    :param timefld: 1d array of the correlation domain
    :returns: 1/e correlation time
    """
    meanvar = mytrapz(data, timefld)
    # The mean needs to be substracted to arrive at
    # the statistical definition of the autocorrelation
    data -= meanvar
    result = np.correlate(data, data, mode="full")
    result = result[int(result.size/2):]
    overe = np.squeeze(np.where(result > np.exp(-1)*result[0]))
    if overe.size <= 1:
        corrtime = 0
    else:
        cort_ind = np.array_split(overe, np.where(np.diff(overe) != 1)[0] + 1)[0][-1]
        corrtime = timefld[cort_ind] - timefld[0]
    return corrtime

## example
##if __name__ == "__main__":
    
