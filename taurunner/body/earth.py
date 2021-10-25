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

def construct_earth(layers: list=[]) -> Body:
    r'''
    Function for making the PREM Earth

    Params
    ______
    layers : Optional list of tuples with radii [km] and densities [gr/cm^3] to add as constant
             density layers on top of the PREM model. for example a 4km ocean layer is [(4, 1.0)]

    Returns
    _______
    earth : TauRunner Earth object
    '''
    r_tot            =  6368.
    layer_boundaries = [1221.5, 3480, 5701, 5771, 5971, 6151, 6346.6, 6356, 6368]
    pparams          = [_ for _ in prem_params] # hacky fix sorry

    for layer in layers:
        r, density = layer
        r_tot += r
        pparams.append((density, 0., 0., 0.))
        layer_boundaries.append(r_tot)
    pparams = list(pparams)
    layer_boundaries = np.array(layer_boundaries, dtype=float) / r_tot
    earth_densities  = [(helper(param), layer) for param, layer in zip(pparams, layer_boundaries)]
    earth            = Body(earth_densities, r_tot, name='PREM_earth')
    return earth
