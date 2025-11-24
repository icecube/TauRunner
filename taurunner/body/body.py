import numpy as np

from abc import ABC, abstractmethod
from scipy.integrate import quad
from typing import Optional, Iterable, Union

from taurunner.utils import units, Callable

class Body(ABC):

    def __init__(
        self,
        density: Iterable[float],
        length: float,
        layer_boundaries: Iterable[float]
    ):
        self._density = [units.gr/units.cm**3*Callable(x) for x in density]
        self._layer_boundaries = layer_boundaries
        self._length = length * units.km

    @property
    def length(self):
        return self._length

    @property
    def density(self):
        return self._density

    @property
    def layer_boundaries(self):
        return self._layer_boundaries

    @abstractmethod
    def get_density(self, x: float):
        pass

    @abstractmethod
    def get_average_density(self, x: float):
        pass
