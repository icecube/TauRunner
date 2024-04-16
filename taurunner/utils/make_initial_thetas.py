import numpy as np

from typing import Union, Tuple

def make_initial_thetas(
    nevents: int,
    theta: Union[float, Tuple[float, float]],
    track_type: str='chord'
) -> np.ndarray:
    r'''
    Creates an array of initial nadir angles
    Params
    ------
    nevents: Number of events
    theta: If float, then initial angle in degrees. If
        tuple, then minimum and maximum angles
    track_type: which type of taurunner.track object
    Returns
    -------
    thetas: Nadir angles in radians
    '''

    if track_type=='radial':
        thetas = np.zeros(nevents)

    elif isinstance(theta, tuple):
        if theta[0]>theta[1]:
            theta = theta[1], theta[0]
        if theta[0]<0 or theta[1]>180:
            raise ValueError('Angles must be between 0 and 180')
        costh_min = np.cos(np.radians(theta[0]))
        costh_max = np.cos(np.radians(theta[1]))
        thetas = np.arccos(np.random.uniform(costh_min, costh_max, nevents))
    else:
        if theta < 0 or theta > 180:
            raise ValueError('Angles must be between 0 and 180')
        thetas = np.radians(np.full(nevents, theta))

    sorter = np.argsort(thetas)
    thetas = thetas[sorter]
    return thetas
