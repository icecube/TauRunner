import numpy as np
from taurunner.utils import is_floatable

def make_initial_thetas(nevents, theta, track_type='chord', rand=None):
    r'''
    Creates an array of initial nadir angles
    Params
    ------
    nevents (int): Number of events
    theta (float or tuple): If float, then initial angle in degrees. If
        tuple, then minimum and maximum angles
    track_type (str, default="chord") : which type of taurunner.track object
    rand (np.random.RandomState): numpy random number generator object
    Returns
    -------
    output (array-like) : Nadir angles in radians
    '''
    if rand is None:
        rand=np.random.RandomState()

    if track_type=='radial':
        thetas = np.zeros(nevents)

    elif is_floatable(theta):
        t = float(theta)
        if t<0 or t>180:
            raise ValueError('Angles must be between 0 and 180')
        thetas = np.radians(np.ones(nevents)*t)
    else:
        if type(theta)==tuple:
            if theta[0]>theta[1]:
                theta = theta[1], theta[0]
            if theta[0]<0 or theta[1]>180:
                raise ValueError('Angles must be between 0 and 180')
            costh_min = np.cos(np.radians(theta[0]))
            costh_max = np.cos(np.radians(theta[1]))
            thetas = np.arccos(rand.uniform(costh_min, costh_max, nevents))
        else:
            raise ValueError('theta sampling %s not suppoorted' % theta)

    return thetas
