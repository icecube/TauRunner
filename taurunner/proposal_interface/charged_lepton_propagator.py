import numpy as np
from proposal.particle import Secondaries

from abc import ABC, abstractmethod

#from .new_proposal_lepton_propagator import make_propagator
#from ..utils import units
from ..particle import Particle
from ..body import Body
from ..track import Track
from ..cross_sections import CrossSections

class ChargedLeptonPropagator(ABC):
    """Class for propagating charged leptons using PROPOSAL"""
    def __init__(self, body: Body):
        self._body = body
        self._propagators = {}

    @property
    def body(self):
        return self._body

    @abstractmethod
    def propagate(
        self,
        particle: Particle,
        track: Track
    )->Secondaries:
        """
        Propagates a charged lepton through a body

        params
        ______
        particle: TauRunner particle to propagate
        body: Body in which to propagate
        track: trajectory along which to propagate

        returns
        _______
        sec: PROPOSAL secodaries
        """
        pass
