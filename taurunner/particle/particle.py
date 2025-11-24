import numpy as np
import proposal as pp
import sys

from typing import Union

from ..utils import units
from ..body import Body
from ..track import Track
from ..cross_sections import CrossSections
from .utils import *


ISOSCALAR_MASS = ((0.9382720813+0.9395654133)/2.) * units.GeV
EMIN = 1e9 # minimum energy allowed in the splines

PROPOSAL_PROPAGATORS = {}

class Particle:
    r'''
    This is the class that contains all relevant 
    particle information stored in an object.
    '''
    def __init__(
        self, 
        ID: int, 
        energy: float, 
        position: float, 
        xs: CrossSections, 
        secondaries: bool = True, 
        no_losses: bool = False
    ):
        r'''

        Params
        ------
        ID: PDG particle identifier
        energy: Initial energy of the particle [eV]
        position: Affine paramter describing the distance along the track of the particle (0<x<1)
        xs: TauRunner CrossSections object
        secondaries: Boolean telling whether to include secondary (mu and e) neutrinos from tau decay
        no_losses: Boolean to turn off charged lepton losses
        '''
        #Set Initial Values
        self.ID = ID
        self.initial_energy = energy
        self.energy = energy
        self.position = position
        #self.chargedposition = 0.0
        self.SetParticleProperties()
        self.secondaries = secondaries
        self.survived = True
        self.basket = []
        self.nCC = 0.0
        self.nNC = 0.0
        self.ntdecay = 0.0
        self.decay_position = 0.0
        self.xs = xs
        self.xs_model = xs.model
        self.losses = not no_losses

        self.nutype = 'nu'
        if ID < 0:
            self.nutype += 'bar'

    def __str__(self):
        return str(self.ID)

    def __repr__(self):
        return f"ID:{self.ID} E:{self.energy/units.GeV} GeV"
        
    def SetParticleProperties(self):
        r'''
        Sets particle properties, either when initializing or after an interaction.
        '''
        if np.abs(self.ID) in [12, 14, 16]:
            self.mass = 0.0 #this is not true.. and it seems to have caused quite the stir.
            self.lifetime = np.inf #this is unclear
        if np.abs(self.ID) == 15:
            self.mass = 1.776*units.GeV
            self.lifetime = 2.9e-13*units.sec
        if np.abs(self.ID) == 13:
            self.mass = 0.105*units.GeV
            self.livetime = 2.2e-6*units.sec

    def GetProposedDepthStep(self):
        r'''
        Calculates the free-streaming column depth of your neutrino based
        on the cross section, and then samples randomly from a log-uniform
        distribution.

        Parameters
        ------------
        p: float
            random number. the free-streaming column depth is scaled by
            the log of this number
        Returns
        -----------
        DepthStep: float
            Column depth to interaction in natural units
        '''
        #Calculate the inverse of the interaction depths.
        #p = self.rand.random_sample()
        p = np.random.random_sample()
        first_piece = (1./self.GetInteractionDepth(interaction='CC'))
        second_piece = (1./self.GetInteractionDepth(interaction='NC'))
        step = (first_piece + second_piece)

        #return the column depth to interaction - weighted by a random number
        return -np.log(p)/step

    def GetTotalInteractionDepth(self):
        return(1./(1./self.GetInteractionDepth('NC')
               + 1./self.GetInteractionDepth('CC')))
        #x = sum([self.GetInteractionDepth(itype)**-1 for itype in "CC NC".split()]) **-1
        #return x

    def GetInteractionDepth(self, interaction: str, proton_fraction: float=0.5):
        r'''
        Calculates the mean column depth to interaction.

        params
        ______
        interaction: str defining the interaction type (CC or NC).

        returns
        _______
        Interaction depth: mean column depth to interaction in natural units
        '''
        if self.ID not in  [12, 14, 16, -12, -14, -16]:
            raise ValueError("Charged leptons shouldn't be in here")
        if not 0 <= proton_fraction <= 1:
            raise ValueError("Proton fraction must be between 0 and 1.")
        if interaction not in "CC NC".split():
            raise ValueError(f"{interaction} is not valid interaction type")

        return ISOSCALAR_MASS/self.xs.total_cross_section(
            self.energy, 
            self.nutype, 
            interaction,
            proton_fraction=proton_fraction
        )
                               

    def GetInteractionProbability(self,ddepth,interaction):
        return 1.-np.exp(-ddepth/self.GetInteractionDepth(interaction))

    def Decay(self):
        if np.abs(self.ID) in [12, 14, 16]:
            raise ValueError("no neutrino decays.. yet")
        if np.abs(self.ID) == 15:
            if self.secondaries:
                # sample branching ratio of tau leptonic decay
                p0 = np.random.uniform(0,1)
                #p0 = self.rand.uniform(0,1)
                if p0 < .18:
                    # sample energy of secondary antinumu
                    #sample = SampleSecondariesEnergyFraction(self.rand.uniform(0,1), antinumu_cdf)
                    sample = SampleSecondariesEnergyFraction(np.random.uniform(0,1), antinumu_cdf)
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : -np.sign(self.ID)*14, "position" : self.position, "energy" : enu})
                elif p0 > .18 and p0 < .36:
                    # sample energy of secondary antinue
                    #sample = SampleSecondariesEnergyFraction(self.rand.uniform(0,1), antinue_cdf)
                    sample = SampleSecondariesEnergyFraction(np.random.uniform(0,1), antinue_cdf)
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : -np.sign(self.ID)*12,  "position" : self.position, "energy" : enu})
            self.energy = self.energy*np.random.choice(TauDecayFractions, p=TauDecayWeights)
            self.ID = np.sign(self.ID) * 16
            self.SetParticleProperties()
            return
        if np.abs(self.ID) in [11, 13]:
            self.survived=False

    def PropagateChargedLepton(
        self,
        body: Body,
        track: Track,
        charged_lepton_propagator
    ):
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
        total_dist = track.x_to_d(1.0 - self.position) * body.length
        if(np.logical_and(np.abs(self.ID) in [13, 14], total_dist / units.km > 100.)):  #muons farther than 100km will not make it
             return
        sec = charged_lepton_propagator.propagate(self, track)
        # Update praticle state
        #self.energy = sec.final_state().energy * units.MeV
        #current_distance = track.x_to_d(self.position)
        #charged_distance = sec.final_state().propagated_distance * units.cm / body.radius
        #self.position = track.d_to_x((current_distance + charged_distance))
        # not adhering to unit protocol here !!!
        #self.chargedposition = sec.final_state().propagated_distance * units.cm / units.km


    def Interact(self, interaction, body=None, track=None, proton_fraction=0.5):
        if np.abs(self.ID) not in [12, 14, 16]:
            raise ValueError('Particle ID not supported by this function') 
        #Sample energy lost from differential distributions
        NeutrinoInteractionWeights = self.xs.differential_cross_section(
            self.energy,
            NeutrinoDifferentialEnergyFractions,
            self.nutype,
            interaction,
            proton_fraction=proton_fraction
        )
        NeutrinoInteractionWeights = np.divide(
            NeutrinoInteractionWeights, 
            np.sum(NeutrinoInteractionWeights)
        )
        z_choice = np.random.choice(
            NeutrinoDifferentialEnergyFractions,
            p=NeutrinoInteractionWeights
        )
        self.energy = z_choice * (self.energy - EMIN) + EMIN

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
