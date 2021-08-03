import numpy as np
import proposal as pp
from importlib.resources import path

import taurunner
from taurunner.modules import units
from .utils import *
from taurunner.resources import secondaries_splines

# TODO move this into utils
# load secondary cdfs
with path(secondaries_splines, 'antinue_cdf.npy') as p:
        nue_path = str(p)
with path(secondaries_splines, 'antinumu_cdf.npy') as p:
        numu_path = str(p)

antinue_cdf = np.load(nue_path)
antinumu_cdf = np.load(numu_path)
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

proton_mass = ((0.9382720813+0.9395654133)/2.)*units.GeV

#Particle object. 
class Particle(object):
    r'''
    This is the class that contains all relevant 
    particle information stored in an object.
    '''
    def __init__(self, ID, energy, incoming_angle, position, seed,
                   xs, proposal_propagator, proposal_lep, secondaries, no_losses):
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
        seed:           int
            Seed corresponding to the random number generator.
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
        self.rand            = np.random.RandomState(seed=seed)
        self.xs              = xs
        self.xs_model        = xs.model
        self.propagator      = proposal_propagator
        self.losses          = not no_losses
        self.lep             = proposal_lep
        
    def SetParticleProperties(self):
        r'''
        Sets particle properties, either when initializing or after an interaction.
        '''
        if self.ID in [12, 14, 16]:
            self.mass = 0.0          #this is not true.. and it seems to have caused quite the stir.
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
            raise ValueError("no neutrino decays.. yet")
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
            #print('uuuuh???')
            self.survived=False

    def PropagateChargedLepton(self, body, track): #description is wrong that should be fixed
        r'''
        propagate taus/mus through medium
        Parameters
        ----------
        objects: list
            List of CasinoEvents that need to have tau losses sampled stochastically.
        xs: str 
            Cross section model to use for the photohadronic losses
        losses: bool
            This can be set to False to turn off energy losses. In this case, the particle decays at rest.
        Returns
        -------
        objects: list
            List of CasinoEvents after losses are calculated
        '''
        #if self.ID==13:
        #  flavor='mu'
        #  lep = pp.particle.DynamicData(pp.particle.MuMinusDef().particle_type)
        #elif self.ID==15: # pragma: no cover
        #  flavor='tau'
        #  lep = pp.particle.DynamicData(pp.particle.TauMinusDef().particle_type)
        lep = self.lep #[self.ID]
        current_km_dist  = track.x_to_d(self.position)*body.radius/units.km
        total_dist       = track.x_to_d(1.-self.position)*body.radius/units.km
        current_density  = body.get_average_density(track.x_to_r(self.position))
        dist_to_prop     = 1e3*(total_dist - current_km_dist)

        if(not self.losses): #easy
            return
        #elif(self.energy/units.GeV <= 1e6):
        #    self.chargedposition = ((self.energy/units.GeV/1e6)*50.)/1e3
        else: # pragma: no cover
            lep_length  = []
            en_at_decay = []
            #need to add support to propagate without decay here (fixed distance propagation)
            lep.energy     = 1e3*self.energy/units.GeV
            rad = body.radius/units.km*1e5
            phi            = 2.*track.theta
            pos_vec   = pp.Vector3D(rad*np.sin(phi)*(1. - self.position), 0, rad*((1. - track.depth + np.cos(phi))*self.position - np.cos(phi)))
            direction = [-np.sin(phi), 0., np.cos(phi) + (1. - track.depth)]
            norm = np.linalg.norm(direction)
            lep.position   = pos_vec
            lep.direction  = pp.Vector3D(direction[0]/norm, 0., direction[2]/norm)
            sec            = self.propagator.propagate(lep) #, dist_to_prop)
            particles      = sec.particles
            final_vec      = (sec.position[-1] - pos_vec)
            lep_length     = final_vec.magnitude() / 1e5
            decay_products = [p for i,p in zip(range(max(len(particles)-3,0),len(particles)), particles[-3:]) if int(p.type) <= 1000000001]
            en_at_decay    = np.sum([p.energy for p in decay_products])
            self.energy    = en_at_decay*units.GeV/1e3
            self.chargedposition  = lep_length
        return

    def Interact(self, interaction, body=None, track=None): #  dist_to_prop=None, current_density=None):
        if self.ID in [12, 14, 16]:
            #Sample energy lost from differential distributions
            NeutrinoInteractionWeights = self.xs.DifferentialOutGoingLeptonDistribution(
                self.energy/units.GeV,
                self.energy*NeutrinoDifferentialEnergyFractions/units.GeV,
                interaction)
            NeutrinoInteractionWeights = np.divide(NeutrinoInteractionWeights, 
                                                   np.sum(NeutrinoInteractionWeights))
            self.energy = self.energy*self.rand.choice(NeutrinoDifferentialEnergyFractions,
                                                       p=NeutrinoInteractionWeights)

            if interaction == 'CC':
                #make a charged particle
                self.nCC += 1
                if(self.ID==16):
                    self.ID = 15
                #elif(self.ID==14):
                #    self.ID = 13
                elif(self.ID in [12, 14]):
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
