import numpy as np
from taurunner.modules import units, Callable
from scipy.integrate import quad

class Body(object):

    def __init__(self,
                 density, 
                 radius: float, 
                 layer_boundaries: list = None, 
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
        layer_boundaries : List which tells the boundaries of different denstiy regions.
                           Must be provided if the density is given as a list, else will be ignored
        # TODO check if name gets used at all....
        name             : Name to give your object
        
        '''

        self.radius = radius*units.km
        self._name  = name
        # Check if body is segmented
        if hasattr(density, '__iter__'):
            if layer_boundaries is None:
                raise RuntimeError('You must specify layer boundaries')
            elif len(layer_boundaries)!=len(density)+1:
                raise RuntimeError('Density and layer_boundaries must have same length.')
            else:
                self._is_layered      = True
                self._density         = [units.gr/units.cm**3*Callable(obj) for obj in density]
                self.layer_boundaries = layer_boundaries
        else:
            self._is_layered      = False
            self._density         = [units.gr/units.cm**3*Callable(density)]
            self.layer_boundaries = np.array([0.0, 1.0])
        self._average_densities()

    def get_density(self, r: float) -> float:
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
        else:
            layer_index = np.digitize(r, self.layer_boundaries)-1
        current_density = self._density[layer_index]
        return current_density(r)

    # TODO do we need this ?
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
