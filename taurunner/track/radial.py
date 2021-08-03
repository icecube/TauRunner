import sys

from .track import Track
from taurunner.modules import doc_inherit

class RadialTrack(Track):

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
