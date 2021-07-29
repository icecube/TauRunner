import numpy as np

import taurunner
from taurunner.modules import units
from .utils import *

# TODO move this into utils
# load secondary cdfs
xs_path = '/'.join(taurunner.__file__.split('/')[:-1]) + '/cross_sections/secondaries_splines/'
#xs_path = os.path.dirname(os.path.realpath(__file__)) + '/cross_sections/secondaries_splines/'
antinue_cdf = np.load(xs_path + 'antinue_cdf.npy')
antinumu_cdf = np.load(xs_path + 'antinumu_cdf.npy')
bins = list(np.logspace(-5,0,101))[:-1] # bins for the secondary splines

def get_sample(u, cdf):
    spl_cdf = iuvs(bins, cdf)
    # check if random sample u is in the range where the spline is defined
    try:
        return (iuvs(bins,cdf-u).roots())[0]
    except:
        if u <= np.min(spl_cdf(bins)): 
            return 1e-3
        elif u == np.max(spl_cdf(bins)):
            return 1

class Particle(object):
    r'''
    This is the class that contains all relevant 
    particle information stored in an object.
    '''
    def __init__(self, ID, energy, incoming_angle, position, index, 
                  seed, chargedposition, xs, secondaries, water_layer=0):
        r'''
        Class initializer. This function sets all initial conditions based 
        on the particle's incoming angle, energy, ID, and position.

        Parameters
        ----------
        ID:    int
            The particle ID as defined in the PDG MC encoding.
        energy:         float
            Initial energy in eV.
        incoming_angle: float
            The incident angle with respect to nadir in radians.
        position:       float
            Position of the particle along the trajectory in natural units (initial should be 0)
        index:          int
            Unique event ID within each run.
        seed:           int
            Seed corresponding to the random number generator.
        chargedposition:    float
            If particle is charged, this is the distance it propagated.
        '''
        #Set Initial Values
        self.ID              = ID
        self.initial_energy  = energy
        self.energy          = energy
        self.position        = position
        self.chargedposition = chargedposition
        self.SetParticleProperties()
        self.secondaries     = secondaries
        self.survived        = True
        self.basket          = []
        self.nCC             = 0
        self.nNC             = 0
        self.ntdecay         = 0
        self.isCC            = False
        self.index           = index
        self.rand            = np.random.RandomState(seed=seed)
        self.xs = xs
        self.xs_model = xs.model

    def SetParticleProperties(self):
        r'''
        Sets particle properties, either when initializing or after an interaction.
        '''
        if self.ID in [12, 14, 16]:
            self.mass = 0.0          #this is not true
            self.lifetime = np.inf   #this is unknown
        if self.ID == 15:
            self.mass = 1.776*units.GeV
            self.lifetime = 2.9e-13*units.sec
        if self.ID == 13:
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

    def GetInteractionDepth(self, interaction):
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
        if self.ID in [12, 14, 16]:
            return proton_mass/(self.xs.TotalNeutrinoCrossSection(self.energy, interaction = interaction))
        if self.ID == 15:
            raise ValueError("Tau interaction length should never be sampled.")

    def GetInteractionProbability(self,ddepth,interaction):
        return 1.-np.exp(-ddepth/self.GetInteractionDepth(interaction))

    def Decay(self):
        if self.ID in [12, 14, 16]:
            raise ValueError("No, you did not just discover neutrino decays..")
        if self.ID == 15:
            if self.secondaries:
                # sample branching ratio of tau leptonic decay
                p0 = self.rand.uniform(0,1)
                if p0 < .18:
                    # sample energy of secondary antinumu
                    sample = get_sample(self.rand.uniform(0,1), antinumu_cdf)
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : 14, "position" : self.position, "energy" : enu})
                elif p0 > .18 and p0 < .36:
                    # sample energy of secondary antinue
                    sample = get_sample(self.rand.uniform(0,1), antinue_cdf)
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : 12,  "position" : self.position, "energy" : enu})
            self.energy = self.energy*self.rand.choice(TauDecayFractions, p=TauDecayWeights)
            self.ID = 16
            self.SetParticleProperties()
            return
        if self.ID == 13:
            self.survived=False


    def Interact(self, interaction):

        if self.ID in [12, 14, 16]:
            #Sample energy lost from differential distributions
            NeutrinoInteractionWeights = self.xs.DifferentialOutGoingLeptonDistribution(
                self.energy/units.GeV,
                self.energy*NeutrinoDifferentialEnergyFractions/units.GeV,
                interaction)
            NeutrinoInteractionWeights = np.divide(NeutrinoInteractionWeights, np.sum(NeutrinoInteractionWeights))
            self.energy = self.energy*self.rand.choice(
                                                       NeutrinoDifferentialEnergyFractions,
                                                       p=NeutrinoInteractionWeights
                                                      )
            
            if interaction == 'CC':
                #make a charged particle
                self.isCC = True
                if(self.ID==16):
                    self.ID = 15
                elif(self.ID==14):
                    self.ID = 13
                self.SetParticleProperties()
                self.nCC += 1

            elif interaction == 'NC':
                #continue being a neutrino
                self.ID = self.ID
                self.SetParticleProperties()
                self.nNC += 1

            return

        elif np.logical_or(self.ID == 15, self.ID == 13):
            raise ValueError("Im not supposed to be here\ntau/mu interactions don't happen here")
