import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import scipy as sp
import pickle
from scipy.interpolate import interp1d
import time
import subprocess
import nuSQUIDSpy as nsq
import CrossSections
from CrossSections import *

units = nsq.Const()
gr = nsq.GlashowResonanceCrossSection()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
tds = nsq.TauDecaySpectra()

earth_model_radii = np.flip([6371., 6367.1774, 6355.7096, 6346.7902, 5960.7076, 5719.8838, 3479.8402, 1221.9578], axis=0)
#earth_model_densities = np.flip([2.6076200000000003, 2.9090978337728903, 3.4256762971354524, 3.900144943911578,
#    4.989777873466213, 11.239115022426343, 12.98114208759252])
earth_model_densities = np.flip([1.360016, 2.6076200000000003, 2.9090978337728903, 3.4256762971354524, 3.900144943911578,
    4.989777873466213, 11.239115022426343, 12.98114208759252], axis=0)

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
        radii = np.sort(radii)
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
        return dists, sects

def check_taundaries(objects, distances):
    r'''
        Propagates Tau leptons while checking boundary conditions
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
    for i, obj in enumerate(objects):
        Ladv = distances[i]*units.km
        #Handle boundary cases, but don't scale by density anywhere
        has_exited = False
        while (obj.region_distance + Ladv >= obj.GetTotalRegionLength()) and (has_exited==False):
    	    try:
		Ladv -= (obj.GetTotalRegionLength() - obj.region_distance)
		obj.position += obj.GetTotalRegionLength() - obj.region_distance
		obj.GetDensityRatio()
		obj.region_distance = 0.
		obj.region_index += 1
	    except:
		has_exited = True
        if has_exited:
	    obj.position = obj.TotalDistance    
        else:
            obj.region_distance += Ladv
            obj.position += Ladv

    return(objects)

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
    e = [obj.energy/units.GeV for obj in objects]                                   #MMC takes initial energy in GeV 
    mult = [obj.GetCurrentDensity()*(units.cm**3)/units.gr/2.7 for obj in objects]  #convert density back to normal (not natural) units 
    sort = sorted(zip(mult, e))
    sorted_mult = np.asarray(zip(*sort)[0])
    sorted_e    = np.asarray(zip(*sort)[1])
    split = np.append(np.append([-1], np.where(sorted_mult[:-1] != sorted_mult[1:])[0]), len(sorted_mult))    
    for i in range(len(split)-1):
        print("SPLIT LENGTH: {}".format(len(sorted_mult[split[i]+1:split[i+1]+1])))
        multis = sorted_mult[split[i]+1:split[i+1]+1]
        eni = sorted_e[split[i]+1:split[i+1]+1]
        #eni_str = " ".join(str(x) for x in eni)
        #arg = '{} {}'.format(eni_str, multis[0]).replace('\n', '')
        #if len(sorted_mult[split[i]+1:split[i+1]+1]) > 1000:
        #    print "ARG = " + arg
        #    print len(arg.split(' '))
        #print len(arg.encode('utf-8'))
        #process = subprocess.Popen('/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/propagate_taus.sh '+arg, stdout=subprocess.PIPE, shell=True)
        eni_str = [str(x) for x in eni]
        eni_str.append(str(multis[0]))
        eni_str.insert(0, '/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/propagate_taus.sh')
        process = subprocess.check_output(eni_str)
        #for line in process.stdout:
        #print(process)
        #print(process.split('\n')[:-1])
        for line in process.split('\n')[:-1]:
            #print(line)
            final_values.append(float(line.replace('\n','')))
    final_energies = np.asarray(final_values)[::2]
    final_distances = np.abs(np.asarray(final_values)[1::2])/1e3
    objects = np.asarray(zip(*sorted(zip(mult, objects)))[1])
    #print 'ARG: ' + arg + '\n\n\n'
    #print final_values
    #print(final_energies.shape)
    #print(final_distances.shape)
    #print(len(objects))
    for i, obj in enumerate(objects):
	obj.energy = final_energies[i]*units.GeV/1e3

    objects = check_taundaries(objects, final_distances)
    return(objects)

Etau = 100.
zz = np.linspace(0.0,1.0,500)[1:-1]
dNTaudz = lambda z: tds.TauDecayToAll(Etau, Etau*z)
TauDecayWeights = np.array(map(dNTaudz,zz))
TauDecayWeights = TauDecayWeights/np.sum(TauDecayWeights)

yy = np.linspace(0.0,1.0,300)[1:-1]

proton_mass = 0.938*units.GeV


class CasinoEvent(object):
    def __init__(self, particle_id, energy, incoming_angle, position, index, CrossSection, seed):
	#Set Initial Values
        self.particle_id = particle_id
        self.energy = energy
        self.position = position
        self.SetParticleProperties()
        self.history = [""]
        self.nCC = 0
        self.nNC = 0
        self.ntdecay = 0
        self.region_index = 0
        self.region_distance = 0.0
        self.CrossSection = CrossSection	
        self.isCC = False
        self.index = index
       	self.rand = np.random.RandomState(seed=seed)

	#Calculate densities along the chord length
	#and total distance
	region_lengths, regions = GetDistancesPerSection(incoming_angle, earth_model_radii)
        densities = []
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
        if self.particle_id == "tau_neutrino":
            self.mass = 0.0
            self.lifetime = np.inf
        if self.particle_id == "tau":
            self.mass = 1.7*units.GeV
            self.lifetime = 1.*units.sec

    def GetParticleId(self):
        return self.particle_id

    def GetLifetime(self):
        return self.lifetime

    def GetMass(self):
        return self.mass

    def GetBoostFactor(self):
        if self.mass > 0.:
            return self.energy/self.mass
        else:
            return np.inf

    def GetProposedDistanceStep(self, density, p):
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
        if self.particle_id == "tau_neutrino":
            # this should be actually divided by average of proton and neutron mass
            return proton_mass/(self.CrossSection.TotalNeutrinoCrossSection(self.energy, interaction = interaction)*density)
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

    def SampleNeutrinoLoss(self):
	return

    def DecayParticle(self):
        if self.particle_id == "tau_neutrino":
            self.history.append("Neutrino decayed???")
            return
        if self.particle_id == "tau":
            self.energy = self.energy*self.CrossSection.rand.choice(zz, p=TauDecayWeights)
            self.particle_id = "tau_neutrino"
            self.SetParticleProperties()
            self.history.append("Tau decayed")
            return

    def InteractParticle(self, interaction):
        if self.particle_id == "tau_neutrino":
  	    #Sample energy lost
            dNdEle = lambda y: self.CrossSection.DifferentialOutGoingLeptonDistribution(self.energy/units.GeV,self.energy*y/units.GeV,interaction = interaction)
            NeutrinoInteractionWeights = map(dNdEle,yy)
	    NeutrinoInteractionWeights = NeutrinoInteractionWeights/np.sum(NeutrinoInteractionWeights)
	    self.energy = self.energy*self.CrossSection.rand.choice(yy, p=NeutrinoInteractionWeights)
           	   
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
		print("i'm not supposed to be here and you should delete these couple of lines")
                event.InteractParticle(interaction = nsq.NeutrinoCrossSections_Current.NC)
    return event

