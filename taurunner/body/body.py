import numpy as np
from taurunner.utils import units, Callable
from scipy.integrate import quad

class Body(object):

    def __init__(self,
                 density, 
                 radius: float, 
                 proton_fraction=0.5,
                 name:str = None
                ):
        r'''
        Object which defines the physical properties of the propagation medium

        Params
        ______
        density          : Object which defines the density of the object. Either float or function
                           that takes a radius (0<=r<=1) and returns a float, or a list of such objects.
                           [gr/cm^3]
        radius           : Radius of the body [meters]
        name             : Name to give your object
        '''

        self.radius    = radius*units.km
        self.radius_km = radius
        self._name     = name
        # Check if body is segmented
        if hasattr(density, '__iter__'):
            self._density    = [units.gr/units.cm**3*Callable(tup[0]) for tup in density]
            layer_boundaries = [tup[1] for tup in density]
            layer_boundaries.insert(0, 0) # the first layer boundary always has to be 0
            self.layer_boundaries = np.array(layer_boundaries)
        else:
            self._density         = [units.gr/units.cm**3*Callable(density)]
            self.layer_boundaries = np.array([0.0, 1.0])
        if hasattr(proton_fraction, '__iter__'):
            self._pfract      = [Callable(tup[0]) for tup in proton_fraction]
            player_boundaries = [tup[1] for tup in proton_fraction]
            player_boundaries.insert(0, 0) # the first layer boundary always has to be 0
            self.player_boundaries = np.array(player_boundaries)
        else:
            self._pfract           = [Callable(proton_fraction)]
            self.player_boundaries = np.array([0.0, 1.0])
        self._average_densities()

    def get_density(self, r: float, right: bool=False) -> float:
        r'''
        Function to get density at an input radius (0<=r<=1)

        Params
        ______
        r : Radius in units of radius of the body, i.e. center is r=0 and edge is r=1

        Returns
        _______
        current_density : Density of the body at the input radius in natural units [eV^4]
        '''
        if r==1:
            layer_index = -1
        elif r==0:
            layer_index = 0
	
        else:
            layer_index = np.digitize(r, self.layer_boundaries, right=right)-1
        current_density = self._density[layer_index]
        return current_density(r)

    def get_average_density(self, r: float) -> float:
        r'''
        Function to return the average density of the layer in which the input radius is.
        Useful for speeding up computation.

        Params
        ______
        r : Radius in units of radius of the body, i.e. center is r=0 and edge is r=1

        Returns
        _______
        avg_dens : Average density in the given layer in natural units [eV^4]
        '''
        if r==1:
            layer_index = -1
        else:
            layer_index = np.digitize(r, self.layer_boundaries)-1
        avg_dens =  self._average_density[layer_index]
        return avg_dens

    def _average_densities(self):
        average_density = []
        for xi,xf in zip(self.layer_boundaries[1:], self.layer_boundaries[:-1]):
            I = quad(self.get_density, xi, xf, full_output=1)
            average_density.append(I[0]/(xf-xi))

        self._average_density = average_density

    def get_proton_fraction(self, r: float, right: bool=False) -> float:
        r'''
        Function to get density at an input radius (0<=r<=1)

        Params
        ______
        r : Radius in units of radius of the body, i.e. center is r=0 and edge is r=1

        Returns
        _______
        current_pfraction : Proton fraction of the body at the input radius in natural units [eV^4]
        '''
        if r==1:
            layer_index = -1
        elif r==0:
            layer_index = 0
        else:
            layer_index = np.digitize(r, self.player_boundaries, right=right)-1
        current_density = self._pfract[layer_index]
        return current_density(r)
