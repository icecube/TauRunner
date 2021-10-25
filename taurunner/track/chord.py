import numpy as np
import sys

from .track import Track
from taurunner.utils import doc_inherit
import proposal as pp

class Chord(Track):

    def __init__(self,
                 depth: float = 0.0,
                 theta: float = 0.0
                ):
        r'''
        depth : Optional argument to specify the depth of the detector
        theta : Incident angle of the particle [radians]
                0 is through the core and \pi/2-\epsilon is skimming
        '''
        Track.__init__(self, depth=depth)
        self.theta = theta
        self.desc  = 'chord'
        self._c    = np.cos(theta)
        self._s    = np.sin(theta)
        self._l1   = np.sqrt(self._c**2+2*depth*self._s**2-depth**2*self._s**2)
        self._l2   = (1-depth)*self._c
        self._t    = self._l1 + self._l2

    def __str__(self):
        desc = (self.theta, self._m)
        return 'theta = %f radians\nm     = %f' % desc

    @doc_inherit
    def d_to_x(self, d: float):
        return d/(self._l1+self._l2)

    @doc_inherit
    def x_to_d_prime(self, x: float):
        return self._l1+self._l2

    @doc_inherit
    def x_to_d(self, x: float):
        return (self._l1+self._l2)*x

    @doc_inherit
    def x_to_cartesian_direction(x: float):
        direction = np.array(-self._s, 0., self._c)
        direction /= np.norm(direction)
        return tuple(direction)

    @doc_inherit
    def r_to_x(self, r: float):
        val1 = (self._l1 - np.sqrt(r**2 - 1 +self._l1**2))/self._t
        val2 = (self._l1 + np.sqrt(r**2 - 1 +self._l1**2))/self._t
        if val2>1 or val1==val2:
            return val1
        else:
            return (val1, val2)

    @doc_inherit
    def x_to_r(self, x: float):
        return np.sqrt(1-2*x*self._l1*(self._l1+self._l2)+x**2*(self._l1+self._l2)**2)
    
    def x_to_r_prime(self, x: float):
        r   = self.x_to_r(x)
        num = -self._l1*self._t + x*self._t**2
        if num==0 and r==0:
            return 0
        else:
            return np.divide(num, r)

    def x_to_pp_dir(self, x: float):
        direction = [-self._s,
                     0, 
                     self._c
                    ]
        dir_vec   = pp.Vector3D(direction[0], 0., direction[2])
        return dir_vec

    def x_to_pp_pos(self, x:float, rad: float):
        #compute direction and position in proposal body
        pos_vec   = pp.Vector3D(self._s*rad*self._t * (1 - x), 
                                0, 
                                -self._c*rad*self._t*(1-x)+(1-self.depth)*rad
                               )
        return pos_vec
