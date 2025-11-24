import numpy as np

from scipy.integrate import quad
from typing import Optional, Iterable, Union

from .body import Body
from taurunner.utils import units, Callable

class SphericalBody(Body):

    def __init__(
        self,
        density: Union[float, Iterable[float]],
        radius: float,
        layer_boundaries: Optional[Iterable[float]]=None,
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

        if hasattr(density, "__iter__"):
            layer_boundaries = [0] + [x[1] for x in density]
            density = [x[0] for x in density]
        else:
            density = [density]
            layer_boundaries = [0.0, 1.0]
        super(SphericalBody, self).__init__(
            density,
            radius,
            layer_boundaries=layer_boundaries
        )

        self._name = name
        self._average_densities()

    @property
    def radius(self):
        return self._length

    @property
    def radius_km(self):
        return self._length / units.km
       
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
        if not (0<=r<=1):
            raise ValueError("r out of range.")

        layer_index = np.digitize(r, self.layer_boundaries, right=True) - 1
        avg_dens =  self._average_density[layer_index]
        return avg_dens

    def _average_densities(self):
        average_density = []
        for xi,xf in zip(self.layer_boundaries[1:], self.layer_boundaries[:-1]):
            I = quad(self.get_density, xi, xf, full_output=1)
            average_density.append(I[0]/(xf-xi))

        self._average_density = average_density
