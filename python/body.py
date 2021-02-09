import numpy as np
from scipy.integrate import quad
import sys
sys.path.append('./modules/')
from physicsconstants import PhysicsConstants
units = PhysicsConstants()

from callable import Callable

class Body():

    def __init__(self, density, radius, layer_boundaries=None):
        r'''
        params
        ______
        density (something) :
        radius (float) : Radius of the body [meters]

        '''

        self.radius  = radius*units.km
        # Check if body is segmented
        if hasattr(density, '__iter__'):
            if layer_boundaries is None:
                raise RuntimeError('You must specify layer boundaries')
            elif len(layer_boundaries)!=len(density)+1:
                raise RuntimeError('Density and layer_boundaries must have same length.')
            else:
                self.density          = [units.gr/units.cm**3*Callable(obj) for obj in density]
                self.layer_boundaries = layer_boundaries
        else:
            self.density = [units.gr/units.cm**3*Callable(density)]
            self.layer_boundaries = np.array([0.0, 1.0])

    def get_density(self, r):
        layer_index     = np.digitize(r, self.layer_boundaries)-1
        current_density = self.density[layer_index]
        return current_density(r)

    def total_column_depth(self, track, safe_mode=False):
        r'''
        params
        ______
        track (Track) :

        returns
        _______
        column_depth (float) :
        '''
        integrand = lambda x: self.get_density(track.x_to_r(x))*track.x_to_r_prime(x)*self.radius
        xx        = np.append(sorted(track.r_to_x(self.layer_boundaries)), 1)
        xx        = xx[np.where(~np.isnan(xx))[0]] # non-nanize
        II        = []
        for xi, xf in zip(xx[:-1], xx[1:]):
            I = quad(integrand, xi, xf)
            if safe_mode:
               if I[1]/I[0] > 1e-4:
                    raise RuntimeError('Error too large')
            II.append(np.abs(I[0]))
        return np.sum(II)/units.gr*units.cm**2
