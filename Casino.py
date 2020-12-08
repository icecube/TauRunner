import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import pythia8
import scipy as sp
import pickle
from scipy.interpolate import interp1d
import time
import subprocess
import nuSQUIDSpy as nsq
from earth import *
import earth

pythia = pythia8.Pythia()

pythia.readString("ProcessLevel:all = off")
pythia.readString("Random:setSeed = on")
pythia.readString("Random:seed = 0")
pythia.init()

units = nsq.Const()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
info = sys.version_info
pyv  = int(info.major)

#cross section tables
######################################
cross_section_path = sys.path[-1]+'../cross_sections/'

if(pyv==3):
    f_NC = np.load(cross_section_path+'NC_table_py3.npy', allow_pickle=True).item()
    f_CC = np.load(cross_section_path+'CC_table_py3.npy', allow_pickle=True).item()

    dsdy_spline_CC = np.load(cross_section_path + 'dsigma_dy_CC_py3.npy', allow_pickle=True).item()
    dsdy_spline_CC_lowe = np.load(cross_section_path + 'dsigma_dy_CC_lowE_py3.npy', allow_pickle=True).item()
    dsdy_spline_NC = np.load(cross_section_path + 'dsigma_dy_NC_py3.npy', allow_pickle=True).item()
    dsdy_spline_NC_lowe = np.load(cross_section_path + 'dsigma_dy_NC_lowE_py3.npy', allow_pickle=True).item()
else:
    f_NC = np.load(cross_section_path+'NC_table.npy', allow_pickle=True).item()
    f_CC = np.load(cross_section_path+'CC_table.npy', allow_pickle=True).item()

    dsdy_spline_CC = np.load(cross_section_path + 'dsigma_dy_CC.npy', allow_pickle=True).item()
    dsdy_spline_CC_lowe = np.load(cross_section_path + 'dsigma_dy_CC_lowE.npy', allow_pickle=True).item()
    dsdy_spline_NC = np.load(cross_section_path + 'dsigma_dy_NC.npy', allow_pickle=True).item()
    dsdy_spline_NC_lowe = np.load(cross_section_path + 'dsigma_dy_NC_lowE.npy', allow_pickle=True).item()

#####################################
#Cross section functions

def TotalNeutrinoCrossSection(enu, xs_model,
                      flavor = nsq.NeutrinoCrossSections_NeutrinoFlavor.tau,
                      neutype = nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
                      interaction = nsq.NeutrinoCrossSections_Current.NC):
        r'''
        Calculates total neutrino cross section. returns the value of sigma_CC (or NC)
        in natural units.
        Parameters
        ----------
        enu:         float
            neutrino energy in eV
        flavor:      nusquids obj
            nusquids object defining neutrino flavor. default is tau
        interaction: nusquids obj
            nusquids object defining the interaction type (CC or NC). default is NC
        Returns
        -------
        TotalCrossSection: float
            Total neutrino cross section at the given energy in natural units.
        '''
        if(xs_model == 'dipole'):
            if(np.log10(enu) < 0.):
                print("Going below a GeV. this region is not supported. godspeed!")
                return 0.
            if(interaction == nsq.NeutrinoCrossSections_Current.NC):
                return((10**f_NC(np.log10(enu/1e9)))*(units.cm)**2)
            else:
                return((10**f_CC(np.log10(enu/1e9)))*(units.cm)**2)
        elif(xs_model == 'CSMS'):
            flavor = nsq.NeutrinoCrossSections_NeutrinoFlavor.tau
            return dis.TotalCrossSection(enu,flavor,neutype,interaction)*(units.cm)**2
        else:
            raise ValueError('Cross section model is not supported. Choose "dipole" or "CSMS".')

def DifferentialOutGoingLeptonDistribution(ein, eout, interaction, xs):
        r'''
        Calculates Differential neutrino cross section. returns the value of d$\sigma$/dy
        in natural units.
        Parameters
        ----------
        ein:         float
            incoming lepton energy in GeV
        eout:         float
            outgoing lepton energy in GeV
        interaction: nusquids obj
            nusquids object defining the interaction type (CC or NC).
        Returns
        -------
        diff: float
            d$\sigma$/d(1-y) at the given energies in natural units where y is the bjorken-y.
        '''
        if(xs=='dipole'):
            if(np.log10(ein) < 0):
                diff = 0
                return diff
            if(interaction==nsq.NeutrinoCrossSections_Current.CC):
                if (np.log10(eout) >= np.log10(ein)):
                    diff = 0.
                    return diff
                elif(np.log10(eout) < 3.):
                    diff = 10**dsdy_spline_CC_lowe(np.log10(ein), np.log10(eout))[0][0]/ein 
                else:
                    diff = 10**dsdy_spline_CC(np.log10(ein), np.log10(eout))[0][0]/ein
            elif(interaction==nsq.NeutrinoCrossSections_Current.NC):
                if (np.log10(eout) >= np.log10(ein)):
                    diff = 0.
                    return diff
                elif( np.log10(eout) <= 3.):
                    diff = 10**dsdy_spline_NC_lowe(np.log10(ein), np.log10(eout))[0][0]/ein
                else:
                    diff = 10**dsdy_spline_NC(np.log10(ein), np.log10(eout))[0][0]/ein
        elif(xs=='CSMS'):
            diff = dis.SingleDifferentialCrossSection(ein*units.GeV, eout*units.GeV, 
                        nsq.NeutrinoCrossSections_NeutrinoFlavor.tau, 
                        nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
                        interaction) 
        return diff

def DoAllCCThings(objects, xs, tau_losses=True):
    r'''
    Calling MMC requires overhead, so handle all MMC calls per iterations
    over the injected events at once
    Parameters
    ----------
    objects: list
        List of CasinoEvents that need to have tau losses sampled stochastically.
    xs: str 
        Cross section model to use for the photohadronic losses
    tau_losses: bool
        This can be set to False to turn off tau energy losses. In this case, the tau decays at rest.
    Returns
    -------
    objects: list
        List of CasinoEvents after Tau Losses have been calculated
    '''
    final_values= []
    efinal, distance = [], []
    e =             [obj[0]/units.GeV for obj in objects]              #MMC takes initial energy in GeV 
    dists =      [1e3*(obj[-2] - obj[1])/units.km for obj in objects]  #distance to propagate in m 
    mult  = [obj[-1]*(units.cm**3)/units.gr/2.7 for obj in objects]    #convert density back to normal (not natural) units
                                                                       #factor of 2.7 because we're scaling the default rho in MMC
    sort = sorted(list(zip(mult, e, dists, objects)))
    sorted_mult = np.asarray(list(zip(*sort))[0])
    sorted_e    = np.asarray(list(zip(*sort))[1])
    sorted_dists = np.asarray(list(zip(*sort))[2])
    sorted_obj = np.asarray(list(zip(*sort))[3])

    if(not tau_losses):
        final_energies = sorted_e
        final_distances = np.zeros(len(sorted_e))
        for i, obj in enumerate(sorted_obj):
            obj[0] = final_energies[i]*units.GeV
            obj[4] = final_distances[i]
        return(sorted_obj)

    split = np.append(np.append([-1], np.where(sorted_mult[:-1] != sorted_mult[1:])[0]), len(sorted_mult))
    
    tau_propagate_path = sys.path[-1]
    if(xs=='dipole'):
        tau_propagate_path+='propagate_taus.sh'
    elif(xs=='CSMS'):
        tau_propagate_path+='propagate_taus_ALLM.sh'
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
            eni_str[kk].insert(0, tau_propagate_path)
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
zz = np.linspace(0.0,1.0,500)[1:-1]
dNTaudz = lambda z: TauDecayToAll(Etau, Etau*z, 0.)
TauDecayWeights = np.array(list(map(dNTaudz,zz)))
TauDecayWeights = np.divide(TauDecayWeights, np.sum(TauDecayWeights))
yy = np.linspace(0.0,1.0,300)[1:-1]

##########################################################
##########################################################
proton_mass = ((0.9382720813+0.9395654133)/2.)*units.GeV

# Event object. where all relevant information about propagation is stored
class CasinoEvent(object):
    r'''
    This is the class that contains all relevant event information stored in an object.
    '''
    def __init__(self, particle_id, energy, incoming_angle, position, index, seed, tauposition, basket, history, water_layer=0, xs_model='dipole', buff=0., body='earth'):
        r'''
        Class initializer. This function sets all initial conditions based on the particle's incoming angle, energy, ID, and position.

        Parameters
        ----------
        particle_id:    string
            The particle can either be "tau" or "tau_neutrino". That is defined here.
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
        tauposition:    float
            If particle is a tau, this is the distance it propagated.
        buff:           float
            If you don't want to simulate to earth emergence, set this buffer to a positive number and simulation will stop short at buff km.
        '''     
        #Set Initial Values
        self.particle_id = particle_id
        self.energy = energy
        self.eini = energy
        self.position = position
        self.SetParticleProperties()
        self.tauposition = tauposition
        self.nCC = 0
        self.nNC = 0
        self.ntdecay = 0
        self.region_index = 0
        self.region_distance = 0.0
        self.isCC = False
        self.index = index
        self.buff = buff
        self.rand = np.random.RandomState(seed=seed)
        self.water_layer = water_layer
        self.basket = basket
        self.history = history
        self.count = 0
        #self.depth = depth
        self.xs_model = xs_model
        if body=='earth':
            earth_model_radii, earth_model_densities = get_radii_densities(self.water_layer)
            #Calculate densities along the chord length
            #and total distance
            region_lengths, regions = GetDistancesPerSection(incoming_angle, earth_model_radii, self.water_layer)
            densities = []
            
            while self.buff > 0:
                if region_lengths[-1] > self.buff:
                    region_lengths[-1] -= self.buff
                    self.buff = 0
                else:
                    self.buff -= region_lengths[-1]
                    region_lengths = np.delete(region_lengths, -1)
                    del regions[-1]
            for i in range(len(region_lengths)):
                region_lengths[i] = region_lengths[i] * units.km
                densities.append(earth_model_densities[regions[i]] * units.gr/(units.cm**3))
       
        elif body=='sun':
            densities      = [152.9 * units.gr/units.cm**3]
            region_lengths = [6.957e5 * units.km, 2*units.km]
        elif body=='christ':
            a_file = open(".ascii_art.txt")
            lines = a_file.readlines()
            for line in lines:
                print(line.strip('\n'))
            quit()
        else:
            print('Body %s is not valid. Only "earth" and "sun" supported presently' % body)

        self.densities = densities
        self.region_lengths = region_lengths
        self.TotalDistance = np.sum(region_lengths)
        cumsum = np.cumsum(self.region_lengths)
        self.region_index = np.min(np.where(cumsum > self.position))
        if(self.position > self.region_lengths[0]):
            self.region_distance = self.position - cumsum[self.region_index -1]
        else:
            self.region_distance = self.position

    def SetParticleProperties(self):
        r'''
        Sets particle properties, either when initializing or after an interaction.
        '''
        if "neutrino" in self.particle_id:
            self.mass = 0.0
            self.lifetime = np.inf
        if self.particle_id == "tau":
            self.mass = 1.7*units.GeV
            self.lifetime = 1.*units.sec

    def GetParticleId(self):
        r'''
        Returns the current particle ID        
        '''
        return self.particle_id

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

    def GetBoostFactor(self):
        r'''
        Returns the current particle's boost factor
        '''
        if self.mass > 0.:
            return self.energy/self.mass
        else:
            return np.inf

    def GetProposedDistanceStep(self, density, p):
        r'''
        Calculates the free-streaming distance of your neutrino based on the density of the medium, and then samples randomly from a log-uniform distribution.

        Parameters
        ------------
        density: float
            medium density in natural units
        p:       float
            random number. the free-streaming distance is scaled by the log of this number
        Returns
        -----------
        DistanceStep: float
            Distance to interaction in natural units
        '''
        #Calculate the inverse of the interaction lengths.
        first_piece = (1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.CC))
        second_piece = (1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.NC))
        step = (first_piece + second_piece)

        #return the distance to interaction - weighted by a random number
        return -np.log(p)/step

    def GetCurrentDensity(self):
        return self.densities[self.region_index]

    def GetTotalRegionLength(self):
        return self.region_lengths[self.region_index]

    def GetTotalInteractionLength(self, density):
        return(1./(1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.NC)
               + 1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.CC)))

    def GetDecayProbability(self,dL):
        boost_factor = self.GetBoostFactor()
        return dL/(boost_factor*self.lifetime)

    def GetInteractionLength(self,density,interaction):
        r'''
        Calculates the mean interaction length.
        
        Parameters
        -----------
        density: float
            medium density in natural units
        interaction: nusquids obj
            nusquids object defining the interaction type (CC or NC).
        Returns
        ----------
        Interaction length: float
            mean interaction length in natural units
        '''
        if "neutrino" in self.particle_id:
            return proton_mass/(TotalNeutrinoCrossSection(self.energy, interaction = interaction, xs_model=self.xs_model)*density)
        if self.particle_id == "tau":
            raise ValueError("Tau interaction length should never be sampled.")

    def GetInteractionProbability(self,dL,density,interaction):
        return 1.-np.exp(-dL/self.GetInteractionLength(density,interaction))

    def GetDensityRatio(self):
        current_density = self.GetCurrentDensity()
        next_density = self.densities[self.region_index + 1]
        return current_density / next_density

    def DecayParticle(self):
        if "neutrino" in self.particle_id:
            #self.history.append("Neutrino decayed???")
            raise ValueError("Dead end.")
        elif self.particle_id == "tau":
            p4vec = pythia8.Vec4(self.eini,0.,0.,self.energy)
            pythia.event.reset()
            pythia.event.append(15,91,0,0,p4vec,self.mass)
            pythia.forceHadronLevel()
            for i in range(pythia.event.size()):
                if not pythia.event[i].isFinal():
                    # if not done decaying, continue
                    continue
                if pythia.event[i].id() == 12:
                    self.basket.append({"id" : "electron_neutrino", "position" : self.position, "energy" : pythia.event[i].e()})
                    # self.history.append("nu_e")
                elif pythia.event[i].id() == 14:
                    self.basket.append({"id" : "muon_neutrino", "position" : self.position, "energy" : pythia.event[i].e()})
                    # self.history.append("nu_mu")
            self.energy = self.energy*self.rand.choice(zz, p=TauDecayWeights)
            self.particle_id = "tau_neutrino"
            self.SetParticleProperties()
            # self.history.append("Tau decayed")
            return

    def check_taundaries(self, tau_position):
        r'''
        Propagates Tau leptons while checking boundary conditions i.e if tau exits earth
        Parameters
        ----------
        objects: list
            list of CasinoEvents
        distances: list
            list of the proposed increased distance step
        Returns
        -------
        objects:
            list of CasinoEvents after propagation
        '''
        if tau_position is None:
            raise Exception('Your tau has not moved. Something is wrong')
        Ladv = (tau_position+0.05*tau_position)*units.km
        has_exited = False
        while (self.region_distance + Ladv >= self.GetTotalRegionLength()) and (has_exited==False):
            try:
                Ladv -= (self.GetTotalRegionLength() - self.region_distance)
                self.position += self.GetTotalRegionLength() - self.region_distance
                self.GetDensityRatio()
                self.region_distance = 0.
                self.region_index += 1
            except:
                has_exited = True
        if has_exited:
            self.position = self.TotalDistance
        else:
            self.region_distance += Ladv
            self.position += Ladv

    def InteractParticle(self, interaction):
        if "neutrino" in self.particle_id:
            #Sample energy lost
            dNdEle = lambda y: DifferentialOutGoingLeptonDistribution(self.energy/units.GeV,self.energy*y/units.GeV, interaction, self.xs_model)
            NeutrinoInteractionWeights = list(map(dNdEle,yy))
            NeutrinoInteractionWeights = np.divide(NeutrinoInteractionWeights, np.sum(NeutrinoInteractionWeights))
            self.energy = self.energy*self.rand.choice(yy, p=NeutrinoInteractionWeights)
                   
            if interaction == nsq.NeutrinoCrossSections_Current.CC:
                self.isCC = True
                self.particle_id = "tau"
                self.SetParticleProperties()
                self.nCC += 1

            elif interaction == nsq.NeutrinoCrossSections_Current.NC:
                self.particle_id = "tau_neutrino"
                self.SetParticleProperties()
                self.nNC += 1
            
            return

        elif self.particle_id == "tau":
            print("Im not supposed to be here")
            raise ValueError("tau interactions don't happen here")


    def PrintParticleProperties(self):
        print("id", self.particle_id, \
              "energy ", self.energy/units.GeV, " GeV", \
              "position ", self.position/units.km, " km")

#This is the propagation algorithm. The MCmeat, if you will.
def RollDice(event):

    while(not np.any((event.position >= event.TotalDistance) or (event.energy <= event.GetMass()) or (event.isCC))):
        if "neutrino" in event.particle_id:
            if(event.energy/units.GeV <= 1e3):
                event.position = event.TotalDistance
                continue
            #Determine how far you're going
            p1 = event.rand.random_sample()
            DistanceStep = event.GetProposedDistanceStep(event.GetCurrentDensity(), p1)
            #Handle boundary conditions and rescale according to column density
            has_exited = False
            while (event.region_distance + DistanceStep >= event.GetTotalRegionLength()) and (has_exited==False):
                try:
                    DistanceStep -= (event.GetTotalRegionLength() - event.region_distance)
                    DistanceStep *= event.GetDensityRatio()
                    event.position += event.GetTotalRegionLength() - event.region_distance
                    event.region_distance = 0.
                    event.region_index += 1
                except:
                    has_exited = True
            if has_exited:
                event.position = event.TotalDistance
                continue
            density = event.GetCurrentDensity()

            #now pick an interaction
            p2 = event.rand.random_sample()
            CC_lint = event.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.CC)
            p_int_CC = event.GetTotalInteractionLength(density) / CC_lint

            if(p2 <= p_int_CC):
                if event.particle_id == "electron_neutrino" or event.particle_id == "muon_neutrino":
                    return None # Kill nues and numus after cc interaction
                event.InteractParticle(nsq.NeutrinoCrossSections_Current.CC)
                event.position += DistanceStep
                event.region_distance += DistanceStep
            else:
                event.InteractParticle(nsq.NeutrinoCrossSections_Current.NC)
                event.position += DistanceStep
                event.region_distance += DistanceStep
            if(event.isCC):
                continue
        elif(event.particle_id == 'tau'):
            event.check_taundaries(event.tauposition)
            if(event.position >= event.TotalDistance):
                return event
                continue
            else:
                event.DecayParticle()
    return event

