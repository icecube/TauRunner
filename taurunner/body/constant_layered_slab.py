import numpy as np

from .body import Body

class ConstantLayeredSlab(body):

    def __init__(
        self,
        density: Union[float, Iterable[float]],
        length: float,
        layer_boundaries: Optional[Iterable[float]]=None
    )

    super(ConstantLayeredSlab, self).__init__(density, length, layer_boundaries)

    def get_density(self, x: float):
        idx = np.digitize(x, self.layer_boundaries, right=True)
        return self.density[idx]

    def get_average_density(self, x: float):
        idx = np.digitize(x, self.layer_boundaries, right=True)
        return self.density[idx]
