import numpy as np
import proposal as pp

from taurunner.body import Body
from taurunner import track
from taurunner.units import units
from taurunner.utils import make_track
from taurunner.particle import Particle
from taurunner.track import Track
from taurunner.cross_sections import CrossSections, XSModel
from taurunner.proposal_interface import ChargedLeptonPropagator

NEUTRINO_IDS = [12, 14, 16, -12, -14, -16]
CHARGED_IDS = [11, 13, 15, -11, -13, -15]
LEPTON_IDS = NEUTRINO_IDS + CHARGED_IDS

def run_MC(
    eini: np.ndarray,
    thetas: np.ndarray,
    body: Body,
    xs: CrossSections,
    seed: int = 0,
    no_secondaries: bool = False,
    flavor: int = 16,
    no_losses: bool = False,
    condition=None,
    depth=0.0,
    track_type='chord'
) -> np.ndarray:
    r'''
    Main simulation code. Propagates an ensemble of initial states and returns the output

    params
    ______
    eini: array containing the initial energies of particles to simulate
    thetas: array containing the incoming angles of particles to simulate
    body: taurunner Body object in which to propagate the particles
    xs: taurunner CrossSections object for interactions
    propagator: PROPOSAL propagator object for charged lepton propagation

    returns
    _______
    output: Array containing the output information of the MC. This includes 
        initial and final energies, incident incoming angles, number of CC and NC interactions,
        and particle type (PDG convention)
    '''

    np.random.seed(seed)
    pp.RandomGenerator.get().set_seed(seed)

    clp = ChargedLeptonPropagator(body, xs)
    nevents = len(eini)
    output = []
    particleIDs = np.full(nevents, flavor, dtype=int)
    xs_model = xs.model
    prev_th = thetas[0]
    prev_track = make_track(prev_th)
    secondary_basket = []
    idxx = []
    my_track  = None
    prv_theta = np.nan
    # Run the algorithm

    # All neutrinos are propagated until exiting as tau neutrino or taus.
    # If secondaries are on, then each event has a corresponding secondaries basket
    # which are propagated all at once in the end.
    for i in range(nevents):
        cur_theta = thetas[i]
        cur_e     = eini[i]

        if (cur_theta!=prv_theta and track_type=='chord') or my_track is None: # We need to make a new track
            my_track = getattr(track, track_type)(theta=cur_theta, depth=depth)
        particle = Particle(
            particleIDs[i],
            cur_e,
            0.0 ,
            xs,
            not no_secondaries,
            no_losses
        )

        out = Propagate(particle, my_track, body, clp, condition=condition)

        if (out.survived==False):
            #this muon/electron was absorbed. we record it in the output with outgoing energy 0
            output.append((cur_e, 0., cur_theta, out.nCC, out.nNC, out.ID, i, out.position))
        else:
            #this particle escaped
            output.append((cur_e, float(out.energy), cur_theta, out.nCC, out.nNC, out.ID, i, out.position))
        if not no_secondaries:
            #store secondaries to propagate later
            secondary_basket.append(np.asarray(out.basket))
            #keep track of parent
            idxx = np.hstack([idxx, [i for _ in out.basket]])
        prv_theta = cur_theta
        del out
        del particle
    idxx = np.asarray(idxx).astype(np.int32)
    if not no_secondaries:
        #make muon propagator
        secondary_basket = np.concatenate(secondary_basket)
        for sec, i in zip(secondary_basket, idxx):
            cur_theta = thetas[i]
            if cur_theta!=prv_theta and track_type=='chord': # We need to make a new track
                my_track = getattr(track, 'chord')(theta=cur_theta, depth=depth)
            sec_particle = Particle(
                sec['ID'], 
                sec['energy'],
                sec['position'],
                xs=xs, 
                secondaries=False, 
                no_losses=False
            )
            sec_out = Propagate(sec_particle, my_track, body, clp, condition=condition)
            if(not sec_out.survived):
                output.append((sec_out.initial_energy, 0.0, cur_theta, sec_out.nCC, sec_out.nNC, sec_out.ID, i, sec_out.position))
            else:
                output.append((sec_out.initial_energy, sec_out.energy, cur_theta, 
    	    		                 sec_out.nCC, sec_out.nNC, sec_out.ID, i, sec_out.position))
            prv_theta = cur_theta
            del sec_particle
            del sec_out

    output = np.array(
        output,
        dtype=[
            ('Eini', float),
            ('Eout',float),
            ('Theta', float),
            ('nCC', int),
            ('nNC', int),
            ('PDG_Encoding', int),
            ('event_ID', int),
            ('final_position', float)
        ]
    )
    output['Theta'] = np.degrees(output['Theta'])
    return output
        

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

            proton_fraction = body.get_proton_fraction(track.x_to_r(particle.position))

            particle.Interact(
                interaction,
                body,
                track,
                proton_fraction=proton_fraction
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
