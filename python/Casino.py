import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import scipy as sp
import pickle
from scipy.interpolate import interp1d
import time
import subprocess
import nuSQUIDSpy as nsq
sys.path.append('./modules')
from physicsconstants import PhysicsConstants
from cross_sections import xs
units = PhysicsConstants()
dis = nsq.NeutrinoDISCrossSectionsFromTables()

def DoAllCCThings(objects, xs, flavor, losses=True):
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
    flavor = objects[0][-1]
    e      = [obj[0]/units.GeV for obj in objects]                    #MMC takes initial energy in GeV 
    dists  = [1e3*(obj[6] - obj[2])/units.km for obj in objects]     #distance to propagate is total distance minus the current position in m 
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
            obj[4] = final_distances[i]
        return(sorted_obj)

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
        nbatch = np.ceil(len(eni)/max_arg)
        eni_str = []
        for i in range(int(nbatch)):
            if i != int(nbatch) - 1:
                eni_str.append(['{} {}'.format(eni[i*max_arg + x], din[i*max_arg + x]) for x in range(max_arg)])
            else:
                eni_str.append(['{} {}'.format(eni[i*max_arg + x], din[i*max_arg + x]) for x in range(len(eni)%max_arg)])

        num_args = len(eni)/max_arg
        if len(eni) % max_arg != 0:
            eni_str.append(["{} {}".format(eni[x], din[x]) for x in range(int(max_arg*num_args), len(eni))])
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
        obj[4] = final_distances[i]
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
        g1 = -(2.*z - 1. - RPion)/(1. - RPion)**2

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
    def __init__(self, ID, flavor, energy, incoming_angle, position, index, 
                  seed, chargedposition, water_layer=0, xs_model='dipole'):
        r'''
        Class initializer. This function sets all initial conditions based 
        on the particle's incoming angle, energy, ID, and position.

        Parameters
        ----------
        ID:    string
            The particle can either be "tau" or "tau_neutrino". That is defined here.
        flavor:         int
            lepton flavor (1:e/2:mu/3:tau)
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
        self.ID = ID
        self.energy = energy
        self.flavor = flavor
        self.position = position
        self.chargedposition = chargedposition
        self.SetParticleProperties()
        self.survived = True
        self.nCC = 0
        self.nNC = 0
        self.ntdecay = 0
        self.isCC = False
        self.index = index
        self.rand = np.random.RandomState(seed=seed)
        self.xs_model = xs_model
          
    def SetParticleProperties(self):
        r'''
        Sets particle properties, either when initializing or after an interaction.
        '''
        if self.ID == "neutrino":
            self.mass = 0.0
            self.lifetime = np.inf
        if self.ID == "tau":
            self.mass = 1.7*units.GeV
            self.lifetime = 1.*units.sec
        if self.ID == "mu":
            self.mass = 0.105*units.GeV
            self.livetime = 1.e-6*units.sec
 
    def GetParticleId(self):
        r'''
        Returns the current particle ID        
        '''
        return self.ID

    def PrintParticleProperties(self):
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

    def GetDecayProbability(self,dL):
        r'''
        Returns the particle's decay probability for a 
        '''
        boost_factor = self.GetBoostFactor()
        return dL/(boost_factor*self.lifetime)

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
        if self.ID == "neutrino":
            return proton_mass/(xs.TotalNeutrinoCrossSection(self.energy, interaction = interaction, xs_model=self.xs_model))
        if self.ID == "tau":
            raise ValueError("Tau interaction length should never be sampled.")

    def GetInteractionProbability(self,ddepth,interaction):
        return 1.-np.exp(-ddepth/self.GetInteractionDepth(interaction))

    def Decay(self):
        if self.ID == "neutrino":
            raise ValueError("No, you did not just discover neutrino decays..")
        if self.ID == "tau":
            self.energy = self.energy*self.rand.choice(TauDecayFractions, p=TauDecayWeights)
            self.ID = "neutrino"
            self.SetParticleProperties()
            return
        if self.ID == "mu":
            self.survived=False

    def Interact(self, interaction):

        if self.ID == "neutrino":

            #Sample energy lost from differential distributions
            dNdEle = lambda y: xs.DifferentialOutGoingLeptonDistribution(self.energy/units.GeV,
                                                           self.energy*y/units.GeV, interaction, self.xs_model)
            NeutrinoInteractionWeights = list(map(dNdEle,NeutrinoDifferentialEnergyFractions))
            NeutrinoInteractionWeights = np.divide(NeutrinoInteractionWeights, np.sum(NeutrinoInteractionWeights))
            self.energy = self.energy*self.rand.choice(NeutrinoDifferentialEnergyFractions,
                                                        p=NeutrinoInteractionWeights)

            if interaction == 'CC':
                #make a charged particle
                self.isCC = True
                if(self.flavor == 3):
                    self.ID = "tau"
                elif(self.flavor== 2):
                    self.ID = "mu"
                self.SetParticleProperties()
                self.nCC += 1

            elif interaction == 'NC':
                #continue being a neutrino
                self.ID = "neutrino"
                self.SetParticleProperties()
                self.nNC += 1

            return

        elif np.logical_or(self.ID == "tau", self.ID == "mu"):
            print("Im not supposed to be here")
            raise ValueError("tau/mu interactions don't happen here")

#This is the propagation algorithm. The MCmeat, if you will.
def Propagate(particle, track, body):
    total_column_depth = track.total_column_depth(body)
    #keep iterating until final column depth is reached or a charged lepton is made
    while(not np.any((particle.position >= 1.) or (particle.isCC))):
        if(particle.ID == 'neutrino'):
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
        elif(np.logical_or(particle.ID == 'tau', particle.ID == 'mu')):
            current_distance=track.x_to_d(particle.position)
            charged_distance = particle.chargedposition*units.km/body.radius
            if(track.d_to_x(current_distance+charged_distance) >=1.):
                particle.position=1.
            else:
                current_distance+=charged_distance
                particle.position=track.d_to_x(current_distance)
            if(particle.position >= 1.):
                return particle
                continue
            else:
                particle.Decay()
    return particle
