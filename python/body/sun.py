import numpy as np
from scipy.interpolate import splev, splrep

from body import Body
from physicsconstants import PhysicsConstants
from callable import Callable
units = PhysicsConstants()

def mass_density_from_model(solar_model_file):
    model   = np.genfromtxt(solar_model_file)
    xx      = model[:,0]
    density = model[:,1]
    tck = splrep(xx, np.log(density))
    func = lambda x: np.exp(splev(x, tck))
    return func

def e_density_from_model(solar_model_file):
    model    = np.genfromtxt(solar_model_file)
    xx       = model[:,0]
    density  = model[:,2]
    tck      = splrep(xx, np.log(density))
    func = lambda x: np.exp(splev(x, tck))
    return func

class Sun(Body):

    def __init__(self, density, radius, edensity, layer_boundaries=None, name='Sun'):

        Body.__init__(self, density, radius, layer_boundaries=layer_boundaries, name=name)
        if self._is_layered:
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

HZ_Sun = make_sun('../../solar_models/Hybrid_model_DM_SUN_HZ.txt')  
LZ_Sun = make_sun('../../solar_models/Hybrid_model_DM_SUN_LZ.txt')  
