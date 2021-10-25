import sys
import proposal as pp

from .track import Track
from taurunner.utils import doc_inherit

class Radial(Track):

    def __init__(self, depth=0.0, theta=0.0):
        
        Track.__init__(self, depth=depth)
        self.theta = 0.0
        self.desc  = 'radial'

    @doc_inherit
    def d_to_x(self, d: float):
        return d/(1-self.depth)

    @doc_inherit
    def r_to_x(self, r: float):
        return r/(1-self.depth)

    @doc_inherit
    def x_to_d(self, x: float):
        return (1-self.depth)*x

    @doc_inherit
    def x_to_d_prime(self, x: float):
        return (1-self.depth)

    @doc_inherit
    def x_to_r(self, x: float):
        return x*(1-self.depth)

    def x_to_pp_dir(self, x: float):
        dir_vec   = pp.Vector3D(0, 0., 1.)
        return dir_vec

    def x_to_pp_pos(self, x:float, rad: float):
        #compute direction and position in proposal body
        phi       = 2.*self.theta
        pos_vec   = pp.Vector3D(0,
                                0, 
                                x*rad*(1-self.depth)
                               )
        return pos_vec
