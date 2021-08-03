import numpy as np
from scipy.interpolate import splev, splrep
from importlib.resources import path

from .body import Body
from taurunner.modules import Callable, units
from taurunner.resources import solar_models

def mass_density_from_model(solar_model_file):
    model   = np.genfromtxt(solar_model_file)
    xx      = model[:,0]
    density = model[:,1]
    tck = splrep(xx, np.log(density))
    return lambda x: np.exp(splev(x, tck))

def e_density_from_model(solar_model_file):
    model    = np.genfromtxt(solar_model_file)
    xx       = model[:,0]
    density  = model[:,2]
    tck      = splrep(xx, np.log(density))
    return lambda x: np.exp(splev(x, tck))

class Sun(Body):

    def __init__(self, density, radius, edensity, layer_boundaries=None, name='Sun'):

        Body.__init__(self, density, radius, layer_boundaries=layer_boundaries, name=name)
        if self._is_layered: # pragma: no cover
            self._edensity = [1./units.cm**3*Callable(obj) for obj in density]
        else:
            self._edensity = [1./units.cm**3*Callable(density)]

    def get_edensity(self, r):
        if r==1:
            layer_index = -1
        else:
            layer_index = np.digitize(r, self.layer_boundaries)-1
        current_edensity = self._edensity[layer_index]
        return current_edensity(r)

def make_sun(model_file):
    density  = mass_density_from_model(model_file)
    edensity = e_density_from_model(model_file)
    return Sun(density, 6.963e5, edensity)

with path(solar_models, 'Hybrid_model_DM_SUN_HZ.txt') as p:
    HZ_path = str(p)
with path(solar_models, 'Hybrid_model_DM_SUN_LZ.txt') as p:
    LZ_path = str(p)

HZ_Sun = make_sun(HZ_path)  
LZ_Sun = make_sun(LZ_path)  
