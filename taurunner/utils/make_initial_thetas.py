import numpy as np
from taurunner.utils import is_floatable

def make_initial_thetas(nevents, theta, track_type='chord', rand=None):

    if rand is None:
        rand=np.random.RandomState()

    if track_type=='radial':
        thetas = np.zeros(nevents)

    elif is_floatable(theta):
        t = float(theta)
        if t<0 or t>90:
            raise ValueError('Angles must be between 0 and 90')
        thetas = np.radians(np.ones(nevents)*t)
    else:
        if type(theta)==tuple:
            if theta[0]<0 or theta[1]>90:
                raise ValueError('Angles must be between 0 and 90')
            costh_min = np.cos(np.radians(theta[0]))
            costh_max = np.cos(np.radians(theta[1]))
            thetas = np.arccos(rand.uniform(costh_min, costh_max, nevents))
        else:
            raise ValueError('theta sampling %s not suppoorted' % theta)

    return thetas
