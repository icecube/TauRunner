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

def chunks(lst, n): # pragma: no cover
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def DoAllCCThings(objects, xs, losses=True):
    r'''
    Calling MMC requires overhead, so handle all MMC calls per iterations
    over the injected events at once
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
    final_values= []
    efinal, distance = [], []
    pid = int(objects[0][-1])
    if pid==14: # pragma: no cover
      flavor='mu'
    elif pid==16: # pragma: no cover
      flavor='tau'
    e      = [obj[0]/units.GeV for obj in objects]                    #MMC takes initial energy in GeV 
    dists  = [1e3*(obj[6] - obj[2])/units.km for obj in objects]      #distance to propagate is total distance minus the current position in m 
    mult   = [obj[-2]*(units.cm**3)/units.gr/2.7 for obj in objects]  #convert density back to normal (not natural) units
    sort         = sorted(list(zip(mult, e, dists, objects)))
    sorted_mult  = np.asarray(list(zip(*sort))[0])
    sorted_e     = np.asarray(list(zip(*sort))[1])
    sorted_dists = np.asarray(list(zip(*sort))[2])
    sorted_obj   = np.asarray(list(zip(*sort))[3])

    if(not losses):
        final_energies = sorted_e
        final_distances = np.zeros(len(sorted_e))
        for i, obj in enumerate(sorted_obj):
            obj[0] = final_energies[i]*units.GeV
            obj[5] = final_distances[i]
        return(sorted_obj)
    else: # pragma: no cover
        split = np.append(np.append([-1], np.where(sorted_mult[:-1] != sorted_mult[1:])[0]), len(sorted_mult))
        propagate_path = os.path.dirname(os.path.realpath(__file__))
        if(xs=='dipole'):
            propagate_path+='/propagate_{}s.sh'.format(flavor)
        elif(xs=='CSMS'):
            propagate_path+='/propagate_{}s_ALLM.sh'.format(flavor)
        else:
            raise ValueError("Cross section model error.")
        for i in range(len(split)-1):
            multis = sorted_mult[split[i]+1:split[i+1]+1]
            eni = sorted_e[split[i]+1:split[i+1]+1]
            din = sorted_dists[split[i]+1:split[i+1]+1]
            max_arg = 500
            eni_str = ['{} {}'.format(e, d) for e,d in list(zip(eni, din))]
            eni_str = list(chunks(eni_str, max_arg))
            for kk in range(len(eni_str)):
                eni_str[kk].append(str(multis[0]))
                eni_str[kk].insert(0, propagate_path)
                process = subprocess.check_output(eni_str[kk])
                for line in process.split(b'\n')[:-1]:
                    final_values.append(float(line.replace(b'\n',b'')))

        final_energies = np.asarray(final_values)[::2]
        final_distances = np.abs(np.asarray(final_values)[1::2])/1e3
        for i, obj in enumerate(sorted_obj):
            obj[0] = final_energies[i]*units.GeV/1e3
            obj[5] = final_distances[i]
        return(sorted_obj)

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
    def __init__(self, ID, energy, incoming_angle, position, index, 
                  seed, chargedposition, water_layer=0, xs_model='dipole', 
                  basket=[]):
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
        self.secondaries     = True
        self.survived        = True
        self.basket          = basket
        self.nCC             = 0
        self.nNC             = 0
        self.ntdecay         = 0
        self.isCC            = False
        self.index           = index
        self.rand            = np.random.RandomState(seed=seed)
        if isinstance(xs_model, CrossSections):
            self.xs = xs_model
            self.xs_model = xs_model.model
        else:
            self.xs_model = xs_model
            self.xs = CrossSections(xs_model)

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
            self.energy = self.energy*self.rand.choice(TauDecayFractions, p=TauDecayWeights)
            self.ID = 16
            self.SetParticleProperties()
            if self.secondaries:
                # sample branching ratio of tau leptonic decay
                p0 = self.rand.uniform(0,1)
                bins = list(np.logspace(-5,0,101))[:-1]
                if p0 < .18:
                    xs_path = os.path.dirname(os.path.realpath(__file__)) + '/cross_sections/secondaries_splines/'
                    cdf = np.load(xs_path + 'antinumu_cdf.npy')
                    # sample energy of tau secondary
                    sample = (iuvs(bins,cdf-np.random.uniform(0,1)).roots())[0]
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : 14, "position" : self.position, "energy" : enu})
                elif p0 > .18 and p0 < .36:
                    xs_path = os.path.dirname(os.path.realpath(__file__)) + '/cross_sections/secondaries_splines/'
                    cdf = np.load(xs_path + 'antinue_cdf.npy')
                    # sample energy of tau secondary
                    sample = (iuvs(bins,cdf-np.random.uniform(0,1)).roots())[0]
                    enu = sample*self.energy
                    # add secondary to basket, prepare propagation
                    self.basket.append({"ID" : 12,  "position" : self.position, "energy" : enu})
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

#This is the propagation algorithm. The MCmeat, if you will.
def Propagate(particle, track, body):
    total_column_depth = track.total_column_depth(body)
    #keep iterating until final column depth is reached or a charged lepton is made
    while(not np.any((particle.position >= 1.) or (particle.isCC))):
        if(particle.ID in [12, 14, 16]):
           #Determine how far you're going
            p1 = particle.rand.random_sample()
            DepthStep = particle.GetProposedDepthStep(p1)

            #now pick an interaction
            p2 = particle.rand.random_sample()
            CC_lint = particle.GetInteractionDepth(interaction='CC')
            p_int_CC = particle.GetTotalInteractionDepth() / CC_lint

            CurrentDepth=track.x_to_X(body, particle.position)
            if(p2 <= p_int_CC):
                if(CurrentDepth+DepthStep >= total_column_depth):
                    particle.position=1.
                    return particle
                else:
                    particle.position=track.X_to_x(body, CurrentDepth+DepthStep)
                    particle.Interact('CC')
            else:
                if(CurrentDepth+DepthStep >= total_column_depth):
                    particle.position=1.
                    return particle
                else:
                    particle.position=track.X_to_x(body, CurrentDepth+DepthStep)
                    particle.Interact('NC')
            if(particle.isCC):
                continue
        elif(np.logical_or(particle.ID == 15, particle.ID == 13)):
            current_distance=track.x_to_d(particle.position)
            charged_distance = particle.chargedposition*units.km/body.radius
            if(track.d_to_x(current_distance+charged_distance) >=1.): # pragma: no cover
                particle.position=1.
            else:
                current_distance+=charged_distance
                particle.position=track.d_to_x(current_distance)
            if(particle.position >= 1.): # pragma: no cover
                return particle
                continue
            else:
                particle.Decay()
    return particle
