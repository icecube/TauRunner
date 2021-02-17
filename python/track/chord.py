import numpy as np

from track import Track

class Chord(Track):

    def __init__(self, depth=0.0, theta=0.0):
        Track.__init__(self, depth=depth)
        self.theta = theta
        self._c    = np.cos(theta)
        self._s    = np.sin(theta)
        self._m    = np.cos(theta) + np.sqrt((1-depth)**2-np.sin(theta)**2)

    def __str__(self):
        desc = (self.theta, self._m)
        return 'theta = %f radians\nm     = %f' % desc

    def d_to_x(self, d, body=None):
        return d/self._m

    def x_to_d(self, x, body=None):
        return self._m*x

    def r_to_x(self, r):
        return (self._c - np.sqrt(r**2- self._s**2))/self._m

    def x_to_r(self, x):
        return np.sqrt(self._s**2 + (self._c-self._m*x)**2)
    
    def x_to_r_prime(self, x):
        return np.abs(self._c-self._m*x)*(-self._m)/self.x_to_r(x)
