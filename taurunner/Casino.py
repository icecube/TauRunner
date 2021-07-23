import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import scipy as sp
import pickle
from scipy.interpolate import interp1d
import time
import subprocess
from scipy.interpolate import InterpolatedUnivariateSpline as iuvs
from taurunner.modules import PhysicsConstants
from taurunner.cross_sections import CrossSections
units = PhysicsConstants()
import proposal as pp

TOL  = 0.0
# load secondary cdfs
xs_path = os.path.dirname(os.path.realpath(__file__)) + '/cross_sections/secondaries_splines/'
antinue_cdf = np.load(xs_path + 'antinue_cdf.npy')
antinumu_cdf = np.load(xs_path + 'antinumu_cdf.npy')
bins = list(np.logspace(-5,0,101))[:-1] # bins for the secondary splines

# Auxiliary function for secondary antinu sampling
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

##########################################################
############## Tau Decay Parameterization ################
##########################################################

RPion = 0.07856**2 
RRho = 0.43335**2 
RA1 = 0.70913**2
BrLepton = 0.18
BrPion = 0.12
BrRho = 0.26
BrA1 = 0.13
BrHad = 0.13

def TauDecayToLepton(Etau, Enu, P):
    z = Enu/Etau
    g0 = (5./3.) - 3.*z**2 + (4./3.)*z**3
    g1 = (1./3.) - 3.*z**2 + (8./3.)*z**3
    return(g0+P*g1)

def TauDecayToPion(Etau, Enu, P):
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RPion - z)  > 0.0):
        g0 = 1./(1. - RPion)
        g1 = -(2.*z - 1. + RPion)/(1. - RPion)**2

    return(g0+P*g1)

def TauDecayToRho(Etau, Enu, P):
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RRho - z) > 0.0):
        g0 = 1./(1. - RRho)
        g1 = -((2.*z-1.+RRho)/(1.-RRho))*((1.-2.*RRho)/(1.+2.*RRho))
    return(g0+P*g1)

def TauDecayToA1(Etau, Enu, P):
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RA1 - z) > 0.0):
        g0 = (1./(1.-RA1))
        g1 = -((2.*z-1.+RA1)/(1.-RA1))*((1.-2.*RA1)/(1.+2.*RA1))
    return(g0 + P*g1)

def TauDecayToHadrons(Etau, Enu, P):
    z = Enu/Etau
    g0=0.
    g1=0.
    if((0.3 - z) > 0.):
        g0 = 1./0.3
    return(g0+P*g1)

def TauDecayToAll(Etau, Enu, P):
    decay_spectra = 0
    decay_spectra+=2.0*BrLepton*TauDecayToLepton(Etau, Enu, P)
    decay_spectra+=BrPion*TauDecayToPion(Etau, Enu, P)
    decay_spectra+=BrRho*TauDecayToRho(Etau, Enu, P)
    decay_spectra+=BrA1*TauDecayToA1(Etau, Enu, P)
    decay_spectra+=BrHad*TauDecayToHadrons(Etau, Enu, P)

    return decay_spectra

Etau = 100.
TauDecayFractions = np.linspace(0.0,1.0,500)[1:-1]
dNTaudz = lambda z: TauDecayToAll(Etau, Etau*z, 0.)
TauDecayWeights = np.array(list(map(dNTaudz,TauDecayFractions)))
TauDecayWeights = np.divide(TauDecayWeights, np.sum(TauDecayWeights))
NeutrinoDifferentialEnergyFractions = np.linspace(0.0,1.0,300)[1:-1]

##########################################################
##########################################################
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
            self.survived=False

    def PropagateChargedLepton(self): #description is wrong that should be fixed
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
        lep = self.lep
        
        if(not self.losses): #easy
            return
        elif(self.energy/units.GeV <= 1e6):
            self.chargedposition = ((self.energy/units.GeV/1e6)*50.)/1e3
        else: # pragma: no cover
            lep_length  = []
            en_at_decay = []
            #need to add support to propagate without decay here (fixed distance propagation)
            lep.energy     = 1e3*self.energy/units.GeV
            #print(lep.energy)
            sec            = self.propagator.propagate(lep) #, distance_to_prop)
            particles      = sec.particles
            lep_length     = sec.position[-1].magnitude() / 1e5     #convert to km
            decay_products = [p for i,p in zip(range(max(len(particles)-3,0),len(particles)), particles[-3:]) if int(p.type) <= 1000000001]
            en_at_decay    = np.sum([p.energy for p in decay_products])
            self.energy    = en_at_decay*units.GeV/1e3
            self.chargedposition  = lep_length
        return

    def Interact(self, interaction,  dist_to_prop=None, current_density=None):
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
                #propagate it
                self.PropagateChargedLepton()
            elif interaction == 'NC':
                #continue being a neutrino
                self.ID = self.ID
                self.SetParticleProperties()
                self.nNC += 1
            return

        elif np.logical_or(self.ID == 15, self.ID == 13):
            raise ValueError("Im not supposed to be here\ntau/mu interactions don't happen here")

#This is the propagation algorithm. The MCmeat, if you will.
def Propagate(particle, track, body):
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
                current_km_dist = track.x_to_d(particle.position)*body.radius/units.km
                #current_density=body.get_density(my_track.x_to_r(out.position))
                current_density  = body.get_average_density(track.x_to_r(particle.position))
                dist_to_prop     = 1e3*(total_distance - current_km_dist)
                particle.Interact('CC', dist_to_prop=dist_to_prop, current_density=current_density)
            else:
                particle.Interact('NC')
        elif(np.logical_or(particle.ID == 15, particle.ID == 13)):
            current_distance=track.x_to_d(particle.position)
            charged_distance = particle.chargedposition*units.km/body.radius
            if(track.d_to_x(current_distance+charged_distance) >=1.): # pragma: no cover
                particle.position=1.
            else:
                current_distance+=charged_distance
                particle.position=track.d_to_x(current_distance)
            if(particle.position >= 1-TOL): # pragma: no cover
                return particle
            else:
                particle.Decay()
    return particle
