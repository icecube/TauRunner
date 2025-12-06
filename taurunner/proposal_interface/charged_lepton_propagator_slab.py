import numpy as np
import proposal as pp

from importlib.resources import path

from .charged_lepton_propagator import ChargedLeptonPropagator
from ..utils import units
from ..particle import Particle
from ..body import Body
from ..track import Track
from ..cross_sections import CrossSections

class ChargedLeptonPropagatorSlab(ChargedLeptonPropagator):
    def __init__(self, slab: Body):
        super(ChargedLeptonPropagatorSlab, self).__init__(slab)

    def propagate(
        self,
        particle: Particle,
        track: Track,
        return_sec=False
    ):

        print(43)
        prv = particle.position
        for (end, medium) in zip(self.body.layer_boundaries,self.body.media):
            if end < particle.position:
                continue
            if (particle.ID, medium) not in self._propagators:
                self._propagators[particle.ID, medium] = make_propagator(particle, medium)
            propagator = self._propagators[(particle.ID, medium)]
            prop_length = (end - prv) * self.body.length
            lep = pp.particle.ParticleState(
                pp.particle.Particle_Type(int(particle.ID - np.sign(particle.ID))),
                pp.Cartesian3D(0.0, 0.0, 0.0),
                pp.Cartesian3D(0.0, 0.0, 1.0),
                particle.energy / units.MeV,
                0.0,
                0.0
            )
            sec = propagator.propagate(lep, int(prop_length / units.cm))
            particle.energy = sec.final_state().energy * units.MeV
            if sec.decay_products():
                break
            prv = end
        if return_sec:
            return sec

particle_def_dict = {
    15: pp.particle.TauMinusDef,
    13: pp.particle.MuMinusDef,
    -15: pp.particle.TauPlusDef,
    -13: pp.particle.MuPlusDef,
}
medium_dict = {
    "StandardRock": pp.medium.StandardRock,
    "Air": pp.medium.Air
}

def make_propagator(particle, medium):

    with path('taurunner.resources.proposal_tables', 'tables.txt') as p:
        tables_path = str(p).split('tables.txt')[0]
    pp.InterpolationSettings.tables_path = tables_path

    particle_def = particle_def_dict[particle.ID]()
    medium = medium_dict[medium.name]()
    collection = pp.PropagationUtilityCollection()

    cuts = pp.EnergyCutSettings(
        np.inf,
        1e-3,
        True
    )
    cross = pp.crosssection.make_std_crosssection(
        particle_def,
        medium,
        cuts,
        True
    )
    create_tables = True
    collection.displacement = pp.make_displacement(cross, create_tables)
    collection.interaction = pp.make_interaction(cross, create_tables)
    collection.time = pp.make_time(cross, particle_def, create_tables)
    collection.decay = pp.make_decay(cross, particle_def, create_tables)

    utility = pp.PropagationUtility(collection=collection)
    geometry = pp.geometry.Sphere(pp.Cartesian3D(0,0,0), 1e20)
    density_dist = pp.density_distribution.density_homogeneous(medium.mass_density)
    prop = pp.Propagator(particle_def, [(geometry, utility, density_dist)])
    return prop

