import proposal as pp
import numpy as np

from taurunner.modules import units
from taurunner.particle import Particle
from taurunner.body import Body
from taurunner.track import Track

#This is the propagation algorithm. The MCmeat, if you will.
def Propagate(particle: Particle, track: Track, body: Body) -> Particle:
    r'''
    Simulate particle along a given track in a body including interactions and energy losses
    Params
    ______
    particle : taurunner Particle object you wish to propagate
    track    : taurunner Track object which defines the geometry of the trajectory
    body     : taurunner Body object which defines the medium in which to simulate
    Returns
    _______
    particle : taurunner Particle object after propagation
    '''
    total_column_depth = track.total_column_depth(body)
    total_distance     = track.x_to_d(1.-particle.position)*body.radius/units.km
    #keep iterating until final column depth is reached or a charged lepton is made
    while(not particle.position >= 1.):
        if(particle.ID in [12, 14, 16]):
            #Determine how far you're going
            p1 = particle.rand.random_sample()
            DepthStep = particle.GetProposedDepthStep(p1)
            CurrentDepth=track.x_to_X(body, particle.position)
            if(CurrentDepth+DepthStep >= total_column_depth):
                particle.position=1.
                return particle
            else:
                particle.position=track.X_to_x(body, CurrentDepth+DepthStep)
            #now pick an interaction
            p2 = particle.rand.random_sample()
            CC_lint = particle.GetInteractionDepth(interaction='CC')
            p_int_CC = particle.GetTotalInteractionDepth() / CC_lint
            if(p2 <= p_int_CC):
                #make and propagate charged lepton
                particle.Interact('CC', body, track)
                particle.PropagateChargedLepton(body, track)
            else:
                particle.Interact('NC')
        elif(np.logical_or(particle.ID == 15, particle.ID == 13)):
            current_distance=track.x_to_d(particle.position)*body.radius
            charged_distance = particle.chargedposition*units.km
            if(track.d_to_x((current_distance+charged_distance)/body.radius) >=1.): # pragma: no cover
                particle.position=1.
            else:
                current_distance+=charged_distance
                particle.position=track.d_to_x(current_distance/body.radius)
                particle.Decay()
    return particle
