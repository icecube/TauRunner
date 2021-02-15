import numpy as np
from scipy.interpolate import splev, splrep

import nuSQUIDSpy as nsq
units = nsq.Const()

class Sun:

    def __init__(self, solar_model_file):

        self.is_discretized   = False
        self.solar_model_file = solar_model_file
        self.model            = np.genfromtxt(solar_model_file)
        self.radius           = 6.957e5 * units.km

        xx       = self.model[:,0]
        mdensity = self.model[:,1]
        edensity = self.model[:,2]
        self.max_x = xx[-1]
        self._mdensity_tck = splrep(xx, np.log(mdensity))
        self._edensity_tck = splrep(xx, np.log(edensity))

    def mass_density(self, x):
        r'''
        Returns mass density at radius $x=r/r_{\odot}$

        params
        ______
        x (float) position along solar radius. Must be between 0 and 1

        returns
        _______
        density (float) mass density at input position [gr/cm^3]
        '''
        if x > self.max_x:
            raise RuntimeError('Radius too large. Maximum radius is x=%f' % self.max_x)
        elif x < 0:
            raise RuntimeError('Radius must be great than 0')
        else:
            return np.exp(splev(x, self._mdensity_tck))

    def electron_density(self, x):
        r'''
        Returns electron density at radius $x=r/r_{\odot}$

        params
        ______
        x (float) position along solar radius. Must be between 0 and 1

        returns
        _______
        density (float) electron density at input position [1/cm^3/NA]
        '''
        if x > self.max_x:
            raise RuntimeError('Radius too large. Maximum radius is x=%f' % self.max_x)
        elif x < 0:
            raise RuntimeError('Radius must be great than 0')
        return np.exp(splev(x, self._edensity_tck))
