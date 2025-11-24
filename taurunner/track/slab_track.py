import numpy as np

from taurunner.body import Body

class SlabTrack:

    def __init__(self):
        pass

    def d_to_x(self, d:float) -> float:
        return d

    def r_to_x(self, r):
        return r

    def total_column_depth(self, body: Body):
        cd, prv = 0, 0
        for (density, x) in zip(body.density, body._layer_boundaries):
            cd += density(None) * (x-prv) * body.length
            prv = x
        return cd

    def x_to_cartesian_direction(self, x):
        return np.array([0, 0, 1])

    def x_to_d(self, x):
        return x

    def x_to_d_prime(self, x):
        return 1.0

    def x_to_r(self, x):
        return x

    def x_to_X(self, body, x):
        cd, prv = 0, 0
        for (density, y) in zip(body.density, body._layer_boundaries):
            if y > x:
                cd += (x - prv) * density(None)
                break
            cd += (y - prv) * density(None)
            prv = y
        return cd * body.length

    def X_to_x(self, body, X):
        prv, x = 0, 0
        for (density, y) in zip(body.density, body._layer_boundaries):
            layer_X = (y - prv) * body.length * density(None)
            if layer_X > X:
                x += X / layer_X * (y - prv)
                break
            X -= layer_X
            x += (y - prv)
            prv = y
        return x


