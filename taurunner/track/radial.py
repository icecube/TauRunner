import sys

from .track import Track
from taurunner.modules import doc_inherit

class RadialTrack(Track):

    @doc_inherit
    def d_to_x(self, d):
        return d/(1-self.depth)

    @doc_inherit
    def r_to_x(self, r):
        return r/(1-self.depth)

    @doc_inherit
    def x_to_d(self, x):
        return (1-self.depth)*x

    @doc_inherit
    def x_to_d_prime(self, x):
        return (1-self.depth)

    @doc_inherit
    def x_to_r(self, x):
        return x*(1-self.depth)
