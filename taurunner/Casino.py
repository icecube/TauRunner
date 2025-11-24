import numpy as np
import proposal as pp

from .utils import units
from .particle import Particle
from .body import Body
from .track import Track
from .cross_sections import CrossSections
from .proposal_interface import ChargedLeptonPropagator

NEUTRINO_IDS = [12, 14, 16, -12, -14, -16]
CHARGED_IDS = [11, 13, 15, -11, -13, -15]
LEPTON_IDS = NEUTRINO_IDS + CHARGED_IDS

class TauRunner:

    def __init__(self, track: Track, body: Body, xs: CrossSections):

        self._track = track
        self._body = body
        self._xs = xs
        self._clp = ChargedLeptonPropagator(body, xs)
        

#This is the propagation algorithm. The MCmeat, if you will.
def Propagate(
    particle: Particle,
    track: Track,
    body: Body,
    charged_lepton_propagator: ChargedLeptonPropagator,
    condition = None
) -> Particle:
    r'''
    Simulate particle along a given track in a body including interactions and energy losses
    params
    ______
    particle : taurunner Particle object you wish to propagate
    track    : taurunner Track object which defines the geometry of the trajectory
    body     : taurunner Body object which defines the medium in which to simulate

    returns
    _______
    particle : taurunner Particle object after propagation
    '''
    total_column_depth = track.total_column_depth(body)
    #keep iterating until final column depth is reached or a charged lepton is made
    stopping_condition = lambda particle: (
        particle.position >= 1. or 
        particle.survived==False
    )

    if condition:
        if hasattr(condition, '__call__'):
            stopping_condition = lambda particle: (particle.position >= 1. or  particle.survived==False or condition(particle))
        elif type(condition)==tuple:
            stopping_condition = lambda particle: (particle.position >= 1. or  particle.survived==False or condition[0](particle, *condition[1]))
        else:
            raise TypeError('Not known how to handle condition argument')

    accumulated_depth = 0
    while(not stopping_condition(particle)):
        if particle.ID not in LEPTON_IDS:
            raise ValueError("Whatcha up to ??")
        
        if particle.ID in NEUTRINO_IDS:

            depth_step = particle.GetProposedDepthStep()
            accumulated_depth += depth_step
            if accumulated_depth >= total_column_depth:
                particle.position = 1.0
                return particle
            else:
                particle.position = track.X_to_x(body, accumulated_depth)

            # I feel like this should all be handled inside particle ?
            # But the is probably cuz particle has too much responsibility
            p2 = np.random.random_sample()

            CC_lint = particle.GetInteractionDepth('CC')
            p_int_CC = particle.GetTotalInteractionDepth() / CC_lint
            interaction = "CC" if p2 <= p_int_CC else "NC"


            particle.Interact(
                interaction,
                body,
                track,
            )

            if interaction=="CC":
                particle.PropagateChargedLepton(
                    body,
                    track,
                    charged_lepton_propagator
                )
        else:
            # This catches rounding errors with PROPOSAL propagation distances
            if 1 - particle.position < 1e-5:
                particle.position = 1
                return particle
            particle.decay_position = particle.position
            particle.Decay()

    return particle
