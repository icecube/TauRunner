import numpy as np
from .body import Body

prem_params = [(13.0885,  0.0,    -8.8381,  0.0),
               (12.5815, -1.2638, -3.6426, -5.5281),
               (7.9565,  -6.4761,  5.5283, -3.0807),
               (5.3197,  -1.4836,  0.0,     0.0),
               (11.2494, -8.0298,  0.0,     0.0),
               (7.1089,  -3.8045,  0.0,     0.0),
               (2.6910,   0.6924,  0.0,     0.0),
               (2.9,      0.0,     0.0,     0.0),
               (2.6,      0.0,     0.0,     0.0),
              ]

def prem_density(r, params):
    return np.polynomial.polynomial.polyval(r, params)
def helper(param):
    func = lambda x: prem_density(x, param)
    return func

def create_earth(layer=None, density=None):
    if(layer==None):
        r_earth          = 6368.
        earth_densities  = [helper(param) for param in prem_params]
        #earth_densities  = [lambda x: prem_density(x, param) for param in prem_params]
        layer_boundaries = np.array([0, 1221, 3480, 5701, 5771, 5971, 6151, 6346.6, 6356, 6368], 
                                    dtype=float) / r_earth
    else:
        r_earth = 6368.+layer
        earth_densities = [helper(param) for param in prem_params.append((density, 0.0, 0.0, 0.0))]
        layer_boundaries = np.array([0, 1221, 3480, 5701, 5771, 5971, 6151, 6346.6, 6356, 6368, r_earth],
                                    dtype=float) / r_earth

    earth = Body(earth_densities, r_earth, layer_boundaries=layer_boundaries, name='PREM_earth')
    return earth
