import numpy as np
import proposal as pp
from importlib.resources import path

import taurunner
from taurunner.utils import units
from taurunner.cross_sections import CrossSections
from .utils import *

from proposal import Propagator
from numpy.random import RandomState
import warnings

proton_mass = ((0.9382720813+0.9395654133)/2.)*units.GeV

ID_2_name = {13:'MuMinusDef', -13:'MuPlusDef', 
             15:'TauMinusDef', -15:'TauPlusDef'}
EMIN = 1e9 # minimum energy allowed in the splines
#Particle object. 
class Particle(object):
    r'''
    This is the class that contains all relevant 
    particle information stored in an object.
    '''
    def __init__(self, 
                 ID:                  int, 
                 energy:              float, 
                 position:            float, 
                 rand:                RandomState,
                 xs:                  CrossSections, 
                 proposal_propagator: Propagator,
                 secondaries:         bool, 
                 no_losses:           bool
                ):
        r'''
        Class initializer. This function sets all initial conditions based 
        on the particle's incoming angle, energy, ID, and position.

        Params
        ------
        ID                  : PDG particle identifier
        energy              : Initial energy of the particle [eV]
        position            : Affine paramter describing the distance along the track of the particle (0<x<1)
        rand                : numpy random number generator
        xs                  : TauRunner CrossSections object
        proposal_propagator : PROPOSAL propagator object for moving charged lepton
        secondaries         : Boolean telling whether to include secondary (mu and e) neutrinos from tau decay
        no_losses           : Boolean to turn off charged lepton losses
        '''
        #Set Initial Values
        self.ID              = ID
        self.initial_energy  = energy
        self.energy          = energy
        self.position        = position
        self.chargedposition = 0.0
        self.SetParticleProperties()
        self.secondaries     = secondaries
        self.survived        = True
        self.basket          = []
        self.nCC             = 0
        self.nNC             = 0
        self.ntdecay         = 0
        self.rand            = rand
        self.xs              = xs
        self.xs_model        = xs.model
        self.propagator      = proposal_propagator
        self.losses          = not no_losses
        if ID > 0:
            self.nutype = 'nu'
        elif ID < 0:
            self.nutype = 'nubar'
        
    def SetParticleProperties(self):
        r'''
        Sets particle properties, either when initializing or after an interaction.
        '''
        if np.abs(self.ID) in [12, 14, 16]:
            self.mass = 0.0          #this is not true.. and it seems to have caused quite the stir.
            self.lifetime = np.inf   #this is unclear
        if np.abs(self.ID) == 15:
            self.mass = 1.776*units.GeV
            self.lifetime = 2.9e-13*units.sec
        if np.abs(self.ID) == 13:
            self.mass = 0.105*units.GeV
            self.livetime = 2.2e-6*units.sec

    def GetParticleId(self):
        r'''
        Returns the current particle ID        
        '''
        return self.ID

    def PrintParticleProperties(self): # pragma: no cover
        print("id", self.ID, \
              "energy ", self.energy/units.GeV, " GeV", \
              "position ", self.position/units.km, " km")

    def GetLifetime(self):
        r'''
        Returns the current particle's lifetime
        '''
        return self.lifetime

    def GetMass(self):
        r'''
        Returns the current particle's mass
        '''
        return self.mass

    def GetProposedDepthStep(self, p):
        r'''
        Calculates the free-streaming column depth of your neutrino based
        on the cross section, and then samples randomly from a log-uniform
        distribution.

        Parameters
        ------------
        p:       float
            random number. the free-streaming column depth is scaled by
            the log of this number
        Returns
        -----------
        DepthStep: float
            Column depth to interaction in natural units
        '''
        #Calculate the inverse of the interaction depths.
        first_piece = (1./self.GetInteractionDepth(interaction='CC'))
        second_piece = (1./self.GetInteractionDepth(interaction='NC'))
        step = (first_piece + second_piece)

        #return the column depth to interaction - weighted by a random number
        return -np.log(p)/step

    def GetTotalInteractionDepth(self):
        return(1./(1./self.GetInteractionDepth(interaction='NC')
               + 1./self.GetInteractionDepth(interaction='CC')))

    def GetInteractionDepth(self, interaction, proton_fraction=0.5):
        r'''
        Calculates the mean column depth to interaction.

        Parameters
        -----------
        interaction: str
            str defining the interaction type (CC or NC).
        Returns
        ----------
        Interaction depth: float
            mean column depth to interaction in natural units
        '''
        if np.abs(self.ID) in [12, 14, 16]:
            return proton_mass/(self.xs.total_cross_section(
                                                            self.energy, 
                                                            self.nutype, 
                                                            interaction,
                                                            proton_fraction=proton_fraction
                                                           )
                               )
        if np.abs(self.ID) == 15:
            raise ValueError("Tau interaction length should never be sampled.")

    def GetInteractionProbability(self,ddepth,interaction):
        return 1.-np.exp(-ddepth/self.GetInteractionDepth(interaction))

    def Decay(self):
        if np.abs(self.ID) in [12, 14, 16]:
            raise ValueError("no neutrino decays.. yet")
        if np.abs(self.ID) == 15:
            if self.secondaries:
                # sample branching ratio of tau leptonic decay
                p0 = self.rand.uniform(0,1)
                if p0 < .18:
                    # sample energy of secondary antinumu
                    sample = SampleSecondariesEnergyFraction(self.rand.uniform(0,1), antinumu_cdf)
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : -np.sign(self.ID)*14, "position" : self.position, "energy" : enu})
                elif p0 > .18 and p0 < .36:
                    # sample energy of secondary antinue
                    sample = SampleSecondariesEnergyFraction(self.rand.uniform(0,1), antinue_cdf)
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : -np.sign(self.ID)*12,  "position" : self.position, "energy" : enu})
            self.energy = self.energy*self.rand.choice(TauDecayFractions, p=TauDecayWeights)
            self.ID = np.sign(self.ID)*16
            self.SetParticleProperties()
            return
        if np.abs(self.ID) in [11, 13]:
            self.survived=False

    def PropagateChargedLepton(self, body, track):
        r'''
        Propagate taus/mus with PROPOSAL along 'track' through 'body'
        Parameters
        ----------
        body: str 
            Cross section model to use for the photohadronic losses
        track: bool
            This can be set to False to turn off energy losses. In this case, the particle decays at rest.
        '''
        if(np.logical_or(not self.losses, np.abs(self.ID) in [11, 12])):
            return
        lep              = pp.particle.DynamicData(getattr(pp.particle, ID_2_name[self.ID])().particle_type)
        total_dist       = track.x_to_d(1.-self.position)*body.radius/units.km
        if(np.logical_and(np.abs(self.ID) in [13, 14], total_dist > 100.)):  #muons farther than 100km will not make it
             return

        lep_length  = []
        en_at_decay = []
        lep.energy     = self.energy/units.MeV
        pos_vec        = track.x_to_pp_pos(self.position, body.radius/units.cm) # radius in cm
        dir_vec        = track.x_to_pp_dir(self.position)
        lep.position   = pos_vec
        lep.direction  = dir_vec
        #propagate
        sec            = self.propagator.propagate(lep, total_dist*1e5) #, dist_to_prop)
        particles      = sec.particles
        #update particle info
        final_vec      = (sec.position[-1] - pos_vec)
        lep_length     = final_vec.magnitude() / 1e5
        decay_products = [p for i,p in zip(range(max(len(particles)-3,0),len(particles)), particles[-3:]) if int(p.type) <= 1000000001]
        en_at_decay    = np.sum([p.energy for p in decay_products])
        if(en_at_decay==0):       #particle reached the border before decaying
            self.energy           = particles[-1].parent_particle_energy*units.MeV
            self.chargedposition  = float(np.ceil(lep_length))
        else:                     #particle decayed before reaching the border
            self.energy           = en_at_decay*units.MeV
            self.chargedposition  = lep_length
        return sec

    def Interact(self, interaction, body=None, track=None, proton_fraction=0.5): #  dist_to_prop=None, current_density=None):
        if np.abs(self.ID) in [12, 14, 16]:
            #Sample energy lost from differential distributions
            NeutrinoInteractionWeights = self.xs.differential_cross_section(self.energy,
                                                                            NeutrinoDifferentialEnergyFractions,
                                                                            self.nutype,
                                                                            interaction,
                                                                            proton_fraction=proton_fraction
                                                    )
            NeutrinoInteractionWeights = np.divide(NeutrinoInteractionWeights, 
                                                   np.sum(NeutrinoInteractionWeights))
            z_choice = self.rand.choice(NeutrinoDifferentialEnergyFractions,
                                                       p=NeutrinoInteractionWeights)
            self.energy = z_choice*(self.energy-EMIN)+EMIN

            if interaction == 'CC':
                #make a charged particle
                self.nCC += 1
                self.ID = np.sign(self.ID)*(np.abs(self.ID)-1)
                if(np.abs(self.ID)==11): #electrons have no chance
                    self.survived=False
                    return
                self.SetParticleProperties()
            elif interaction == 'NC':
                #continue being a neutrino
                self.ID = self.ID
                self.SetParticleProperties()
                self.nNC += 1
            
            return
        else:
            raise ValueError('Particle ID not supported by this function') 
