import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import scipy as sp
import pickle
from scipy.interpolate import interp1d
import time
import subprocess
import nuSQUIDSpy as nsq

units = nsq.Const()
gr = nsq.GlashowResonanceCrossSection()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
tds = nsq.TauDecaySpectra()


#cross sections patched with nuSQUIDS 
#temporary until EHE cross sections are added to nuSQUIDS
######################################
cross_section_path = sys.path[-1]+'../cross_sections/'
charged = np.load(cross_section_path + 'nuXS_CC_8-16.npy')
neutral = np.load(cross_section_path + 'nuXS_NC_8-16.npy')
energies = np.logspace(17, 25, 500)
energies = (energies[:-1] + energies[1:])/2

log_e = np.log10(energies)
log_XS_CC = np.log10(charged)
log_XS_NC = np.log10(neutral)

f_NC = interp1d(log_e, log_XS_NC)
f_CC = interp1d(log_e, log_XS_CC)

dsdy_spline_CC = np.load(cross_section_path + 'dsigma_dy_CC.npy').item()
dsdy_spline_CC_lowe = np.load(cross_section_path + 'dsigma_dy_CC_lowE.npy').item()
dsdy_spline_NC = np.load(cross_section_path + 'dsigma_dy_NC.npy').item()
dsdy_spline_NC_lowe = np.load(cross_section_path + 'dsigma_dy_NC_lowE.npy').item()

#####################################

def TotalNeutrinoCrossSection(enu,
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

	if(np.log10(enu) > log_e[0]):

	    if(interaction == nsq.NeutrinoCrossSections_Current.NC):
		return((10**f_NC(np.log10(enu)))*(units.cm)**2)
            else:
                return((10**f_CC(np.log10(enu)))*(units.cm)**2)
	else:
	    return dis.TotalCrossSection(enu,flavor,neutype,interaction)*(units.cm)**2

def DifferentialOutGoingLeptonDistribution(ein, eout, interaction):
	r'''
	Calculates Differential neutrino cross section. returns the value of d$\sigma$/dy
	in natural units.
	Parameters
	----------
	ein:         float
	    incoming lepton energy in eV
	eout:         float
	    outgoing lepton energy in eV
	interaction: nusquids obj
	    nusquids object defining the interaction type (CC or NC).
	Returns
	-------
	diff: float
	    d$\sigma$/dy at the given energies in natural units.
	'''

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

        return diff

#########################################################
#####  Earth Model Handling. Here we define PREM ########
####   and handling of ice layers/detector       ########
#########################################################

earth_model_radii = np.flip([6371., 6369., 6367.1774, 6355.7096, 6346.7902, 5960.7076, 5719.8838, 3479.8402, 1221.9578], axis=0)
earth_model_densities = np.flip([1.0, 1.360016, 2.6076200000000003, 2.9090978337728903, 3.4256762971354524, 3.900144943911578,
    4.989777873466213, 11.239115022426343, 12.98114208759252], axis=0)

def outer_shell(theta, r1, r2):
    r'''
    Manually add an outer layer of ice to the detector.
    This is especially important for earth skimming trajectories
    Parameters
    -----------
    theta: float
        incoming nadir
    r1: float
        radius of inner shell
    r2: float
        radius of outer shell
    Returns
    --------
    l:  float
        distance travelled in the new layer
    '''
    if theta == np.pi / 2.:
        return np.sqrt(r2**2 - r1**2)
    else:
        return -1*np.cos(theta) * (r1 - np.sqrt(r1**2 + (r2**2 - r1**2) / (np.cos(theta)**2))) 

def GetDistancesPerSection(theta, radii):
        r'''
        Calculates distance traveled in each uniform density region,
        given in order with index of density region
        Parameters
        ----------
        radii: list
            radii of earth regions of different uniform density
        Returns
        -------
        dists_by_sect: list
            distance traveled in each section of uniform density,
            along with the index of the region in an additional dimension
        '''
        radii = np.sort(radii)[:-1]
        closest_approach = radii[-1] * np.sin(theta)
        smallest_sect = np.min(np.where(closest_approach < radii))
        dists = []
        for i in range(smallest_sect, len(radii)):
            side_length = np.sqrt(radii[i]**2 - closest_approach**2)
            dists.append(side_length - np.sum(dists))
        dists = np.flip(dists, axis=0)
        dists[-1] = 2 * dists[-1]
        dists = np.append(dists, np.flip(dists[:-1], axis=0))
        sects = list(range(len(radii) - 1, smallest_sect, -1)) + \
                    list(range(smallest_sect, len(radii)))
        dists = np.insert(dists, 0, outer_shell(theta, radii[-1], 6371.))
        sects.insert(0, len(radii) ) #Set outer shell density to that of ice
        return dists, sects

###########################################
##### Earth Model handling ends here ######
###########################################

def DoAllCCThings(objects):
    r'''
    Calling MMC requires overhead, so handle all MMC calls per iterations
    over the injected events at once
    Parameters
    ----------
    objects: list
        List of CasinoEvents that need to have tau losses sampled
    Returns
    -------
    objects: list
        List of CasinoEvents after Tau Losses have been calculated    
    '''
    final_values= []
    efinal, distance = [], []
    e = [obj[0]/units.GeV for obj in objects]                                #MMC takes initial energy in GeV 
    mult = [obj[-1]*(units.cm**3)/units.gr/2.7 for obj in objects]           #convert density back to normal (not natural) units - 
									     # - factor of 2.7 because thats the default rho in MMC - mult will scale it 
    sort = sorted(zip(mult, e, objects))
    sorted_mult = np.asarray(zip(*sort)[0])
    sorted_e    = np.asarray(zip(*sort)[1])
    sorted_obj = np.asarray(zip(*sort)[2])
    split = np.append(np.append([-1], np.where(sorted_mult[:-1] != sorted_mult[1:])[0]), len(sorted_mult))
    for i in range(len(split)-1):
        multis = sorted_mult[split[i]+1:split[i+1]+1]
        eni = sorted_e[split[i]+1:split[i+1]+1]
        max_arg = 5000
        eni_str = [[str(eni[y*max_arg + x]) for x in range(max_arg)] for y in range(len(eni)/max_arg)]
        num_args = len(eni)/max_arg
        if len(eni) % max_arg != 0:
            eni_str.append([str(eni[x]) for x in range(max_arg*num_args, len(eni))])
        for kk in range(len(eni_str)):
            eni_str[kk].append(str(multis[0]))
            eni_str[kk].insert(0, '/data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/propagate_taus.sh')
            process = subprocess.check_output(eni_str[kk])
            for line in process.split('\n')[:-1]:
                final_values.append(float(line.replace('\n','')))
    final_energies = np.asarray(final_values)[::2]
    final_distances = np.abs(np.asarray(final_values)[1::2])/1e3
    for i, obj in enumerate(sorted_obj):
        obj[0] = final_energies[i]*units.GeV/1e3
        obj[4] = final_distances[i]
#    objects = check_taundaries(sorted_obj, final_distances)
    return(sorted_obj)

##########################################################
############## Tau Decay Parameterization ################
##########################################################

Etau = 100.
zz = np.linspace(0.0,1.0,500)[1:-1]
dNTaudz = lambda z: tds.TauDecayToAll(Etau, Etau*z)
TauDecayWeights = np.array(map(dNTaudz,zz))
TauDecayWeights = TauDecayWeights/np.sum(TauDecayWeights)

yy = np.linspace(0.0,1.0,300)[1:-1]

##########################################################
proton_mass = 0.938*units.GeV

class CasinoEvent(object):
     r'''
     This is the class that contains all relevant event information stored in an object.
     '''
    def __init__(self, particle_id, energy, incoming_angle, position, index, seed, tauposition, buff=0.):
        r'''
        Class initializer. This function sets all initial conditions based on the particle's incoming angle, energy, ID, and position.

        Parameters
        ----------
        particle_id:    string
            The particle can either be "tau" or "tau_neutrino". That is defined here.
	energy:         float
	    Initial energy in eV.
	incoming_angle: float
	    The incident angle with respect to nadir in radians
	position:       float
	    Position of the particle along the trajectory in natural units (initial should be 0)
        index:          int
	    This parameter is purely for bookkeeping purposes.
        seed:           int
	    Seed corresponding to the random number generator. to ensure both randomness and reproducibility.
	tauposition:    float
	    I don't remember what this is.
	buff:           float
	    If you don't want to simulate to earth emergence, set this buffer to a positive number and simulation will stop short at buff km.
        '''	
        #Set Initial Values
        self.particle_id = particle_id
        self.energy = energy
        self.position = position
        self.SetParticleProperties()
        self.tauposition = tauposition
#        self.history = [""]
        self.nCC = 0
        self.nNC = 0
        self.ntdecay = 0
        self.region_index = 0
        self.region_distance = 0.0
        self.isCC = False
        self.index = index
        self.buff = buff
       	self.rand = np.random.RandomState(seed=seed)

        #Calculate densities along the chord length
        #and total distance
        region_lengths, regions = GetDistancesPerSection(incoming_angle, earth_model_radii)
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
        self.region_lengths = region_lengths
        self.densities = densities
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
        if self.particle_id == "tau_neutrino":
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
        #Calculate the inverse of the interaction lengths (The lengths should be added reciprocally)
        first_piece = (1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.CC))
        second_piece = (1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.NC))
        step = (first_piece + second_piece)

        #return the distance to interaction - weighted by a sampled random number
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
        if self.particle_id == "tau_neutrino":
            # this should be actually divided by average of proton and neutron mass
            return proton_mass/(TotalNeutrinoCrossSection(self.energy, interaction = interaction)*density)
        if self.particle_id == "tau":
            # here we need the total tau cross section
            #return TotalNeutrinoCrossSection(self.energy)*density/(proton_mass)
            return 0.01*units.cm

    def GetInteractionProbability(self,dL,density,interaction):
        return 1.-np.exp(-dL/self.GetInteractionLength(density,interaction))

    def GetDensityRatio(self):
        current_density = self.GetCurrentDensity()
        next_density = self.densities[self.region_index + 1]
        return current_density / next_density

    def DecayParticle(self):
        if self.particle_id == "tau_neutrino":
            #self.history.append("Neutrino decayed???")
            return
        if self.particle_id == "tau":
            self.energy = self.energy*self.rand.choice(zz, p=TauDecayWeights)
            self.particle_id = "tau_neutrino"
            self.SetParticleProperties()
            #self.history.append("Tau decayed")
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
        Ladv = tau_position*units.km
        #Handle boundary cases, but don't scale by density anywhere
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
        if self.particle_id == "tau_neutrino":
  	    #Sample energy lost
            dNdEle = lambda y: DifferentialOutGoingLeptonDistribution(self.energy/units.GeV,self.energy*y/units.GeV,interaction = interaction)
            NeutrinoInteractionWeights = map(dNdEle,yy)
	    NeutrinoInteractionWeights = NeutrinoInteractionWeights/np.sum(NeutrinoInteractionWeights)
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
            Efin, Ladv = self.SampleFinalTauParams()
            Ladv *= units.km
            self.energy = Efin*units.GeV
            #Handle boundary cases, but don't scale by density anywhere
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
            self.region_distance += Ladv
            self.position += Ladv
            if (self.position >= self.TotalDistance):
                return
            else:
                self.DecayParticle()
                return

    def PrintParticleProperties(self):
        print "id", self.particle_id, \
              "energy ", self.energy/units.GeV, " GeV", \
              "position ", self.position/units.km, " km"

def RollDice(event):

    while(not np.any((event.position >= event.TotalDistance) or (event.energy <= event.GetMass()) or (event.isCC))):
        if(event.particle_id == 'tau_neutrino'):
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

