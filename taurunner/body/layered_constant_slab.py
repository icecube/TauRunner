import numpy as np

from typing import Union, Iterable, Optional

from .body import Body

from enum import Enum

class Media(Enum):
    Air = 1
    StandardRock = 2

class LayeredConstantSlab(Body):

    def __init__(
        self,
        density: Union[float, Iterable[float]],
        media: Union[str, Iterable[str]],
        length: float,
        layer_boundaries: Optional[Iterable[float]]=None
    ):

        if not hasattr(density, "__iter__"):
            if type(media)!=str:
                raise ValueError("Bad media")
            media = [media]
            layer_boundaries = [1.0]
            density = [density]
        elif layer_boundaries is None:
            raise ValueError("Must specify layer boundaries")

        if not (len(density)==len(media)==len(layer_boundaries)):
            raise ValueError("Incompatible layers")
        super(LayeredConstantSlab, self).__init__(density, length, layer_boundaries)
        self._media = [getattr(Media, x) for x in media]

    @property
    def media(self):
        return self._media

    def get_density(self, x: float):
        if not (0 <= x <= 1):
            raise ValueError(f"{x} not in range [0,1]")
        layer_idx = np.digitize(x, self._layer_boundaries, right=True)
        return self._density[layer_idx]

    def get_average_density(self, x: float):
        return self.get_density(x)
