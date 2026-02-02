import numpy as np
import proposal as pp

from .charged_lepton_propagator import ChargedLeptonPropagator
from .new_proposal_lepton_propagator import make_propagator
from ..utils import units
from ..particle import Particle
from ..body import Body
from ..track import Track
from ..cross_sections import CrossSections

class SphericalBodyChargedLeptonPropagator(ChargedLeptonPropagator):
    """Class for propagating charged leptons using PROPOSAL"""
    def __init__(
            self,
            body: Body,
            xs: CrossSections,
        ):
        super(SphericalBodyChargedLeptonPropagator, self).__init__(
            body,
            xs
        )

    def propagate(self, particle: Particle, body: Body, track: Track):
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
        if particle.ID not in self._propagators:
            self._propagators[particle.ID] = make_propagator(
                particle.ID,
                self._body,
                self._xs
            )
        propagator = self._propagators[particle.ID]
        
        total_dist = track.x_to_d(1.0 - particle.position) * body.radius
        lep = pp.particle.ParticleState(
            pp.particle.Particle_Type(int(particle.ID - np.sign(particle.ID))),
            track.x_to_pp_pos(particle.position, body.radius/units.cm),
            track.x_to_pp_dir(particle.position),
            particle.energy / units.MeV,
            0.0,
            0.0
        )

        # `int` here rounds to nearest cm and stops boundary segfaults
        sec = propagator.propagate(lep, int(total_dist / units.cm))
        # I feel like I shouldn't be returning a PROPOSAL object if this is
        # truly supposed to be an interface
        return sec
