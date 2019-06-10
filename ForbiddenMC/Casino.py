#!/usr/bin/env python

import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import scipy as sp
import nuSQUIDSpy as nsq
import pickle
import argparse
from scipy.interpolate import interp1d
import time
import subprocess

#finaltaus = 0
#earththing = 0
#initialization = 0
#starttime = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('-s',dest='seed',type=int,help='just an integer seed to help with output file names')
parser.add_argument('-n', dest='nevents', type=float, help='how many events do you want?')
parser.add_argument('-gzk', dest='gzk', default=False, action='store_true', help='do you want to propagate the GZK flux? if so, raise this flag, and raise your flag, and raise your flag, and raise it.. when i get older...')
parser.add_argument('-e', dest='energy', type=float, help='if you want to simulate a specific energy, pass it here in GeV')
parser.add_argument('-t', dest='theta', type=float, help='zenith angle in radians where 0 is through the core')
args = parser.parse_args()

#base_path = os.getcwd()+'/'
base_path = '/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/'
cross_section_path = '/data/user/isafa/ANITA/features/TauDragon/cross_sections/'

if (not (args.seed or args.nevents)):
    raise RuntimeError('You must specify a seed (-s) and number of events to simulate (-n)') 

seed = args.seed
nevents = int(args.nevents)
isgzk = args.gzk

units = nsq.Const()
gr = nsq.GlashowResonanceCrossSection()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
tds = nsq.TauDecaySpectra()

rand = np.random.RandomState(seed=seed)

if(isgzk):
  # sample initial energies and incoming angles
  #cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
  #thetas = np.arccos(cos_thetas)
  thetas = np.zeros(nevents)

  gzk_cdf = np.load(base_path+'gzk_cdf_phi_spline.npy').item()
  cdf_indices = rand.uniform(size=nevents)
  #eini = gzk_cdf(cdf_indices)*units.GeV
  eini = np.ones(nevents)*1e9*units.GeV
else:
    eini = np.ones(nevents)*args.energy*units.GeV
    thetas = np.ones(nevents)*args.theta


#cross sections from 1e8 - 1e16 GeV patched with nuSQUIDS 
#temporary until EHE cross sections are added to nuSQUIDS

######################################
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
dsdy_spline_NC = np.load(cross_section_path + 'dsigma_dy_NC.npy').item()

#####################################

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

#def DifferentialOutGoingLeptonDistribution(enu_in,enu_out,
#                                       flavor = nsq.NeutrinoCrossSections_NeutrinoFlavor.tau,
#                                       neutype = nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
#                                       interaction = nsq.NeutrinoCrossSections_Current.NC
#                                    ):
#    diff = dis.SingleDifferentialCrossSection(enu_in,enu_out,flavor,neutype,interaction)
#    return diff

#interaction = 'CC'

def DifferentialOutGoingLeptonDistribution(ein, eout, interaction):
    if(interaction==nsq.NeutrinoCrossSections_Current.CC):
        if (np.log10(eout) < 3. or np.log10(eout) > np.log10(ein)):
            diff = 0.
        else:
            diff = 10**dsdy_spline_CC(np.log10(ein), np.log10(eout))[0][0]/ein
    elif(interaction==nsq.NeutrinoCrossSections_Current.NC):
        if (np.log10(eout) < 3. or np.log10(eout) > np.log10(ein)):
            diff = 0.
        else:
            diff = 10**dsdy_spline_NC(np.log10(ein), np.log10(eout))[0][0]/ein
        
    return diff

def check_taundaries(objects, distances):
    for i, obj in enumerate(objects):
        Ladv = distances[i]*units.km
        #Handle boundary cases, but don't scale by density anywhere
        has_exited = False
        while (obj.region_distance + Ladv >= obj.GetRegionDistance()) and (has_exited==False):
    	    try:
		Ladv -= (obj.GetRegionDistance() - obj.region_distance)
		obj.position += obj.GetRegionDistance() - obj.region_distance
		obj.GetDensityRatio()
		obj.region_distance = 0.
		obj.region_index += 1
	    except:
		has_exited = True
        if has_exited:
	    obj.position = obj.TotalDistance
        obj.region_distance += Ladv
        obj.position += Ladv

    return(objects)

def DoAllCCThings(objects):
    final_values= []
    efinal, distance = [], []
    e = [obj.energy/units.GeV for obj in objects]                                   #MMC takes initial energy in GeV 
    mult = [obj.GetCurrentDensity()*(units.cm**3)/units.gr/2.7 for obj in objects]  #convert density back to normal (not natural) units 
    sort = sorted(zip(mult, e))
    sorted_mult = np.asarray(zip(*sort)[0])
    sorted_e    = np.asarray(zip(*sort)[1])
    split = np.append(np.append([-1], np.where(sorted_mult[:-1] != sorted_mult[1:])[0]), len(sorted_mult))    
    print(split) 
    for i in range(len(split)-1):
        multis = sorted_mult[split[i]+1:split[i+1]+1]
        eni = sorted_e[split[i]+1:split[i+1]+1]
        print(len(eni))
	arg = '{} {}'.format(eni, multis[0]).replace('\n', '').replace('[','').replace(']','')
	process = subprocess.Popen('/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/propagate_taus.sh '+arg, stdout=subprocess.PIPE, shell=True)
        for line in process.stdout:
            final_values.append(float(line.replace('\n','')))
    print('im out!!!!')
    final_energies = np.asarray(final_values)[::2]
    final_distances = np.abs(np.asarray(final_values)[1::2])/1e3
    print(len(final_energies)+len(final_distances))
    objects = np.asarray(zip(*sorted(zip(mult, objects)))[1])
    print(final_energies.shape)
    print(final_energies) 
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
    def __init__(self, particle_id, energy, position, distances, densities, index):
        ## need to add tree in an efficient way

        self.particle_id = particle_id
        self.energy = energy
        self.position = position
        self.SetParticleProperties()
        self.history = [""]
        self.distances = distances
        self.densities = densities
        self.nCC = 0
        self.nNC = 0
        self.ntdecay = 0
        self.region_index = 0
        self.region_distance = 0.0
        self.TotalDistance = np.sum(self.distances)
        self.isCC = False
        self.index = index

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

    def GetRegionDistance(self):
        return self.distances[self.region_index]

    def GetTotalInteractionLength(self, density):
        return(1./(1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.NC)
               + 1./self.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.CC)))

    def GetDecayProbability(self,dL):
        boost_factor = self.GetBoostFactor()
        return dL/(boost_factor*self.lifetime)

    def GetInteractionLength(self,density,interaction):
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
            self.history.append("Neutrino decayed???")
            return
        if self.particle_id == "tau":
            self.energy = self.energy*rand.choice(zz, p=TauDecayWeights)
            self.particle_id = "tau_neutrino"
            self.SetParticleProperties()
            self.history.append("Tau decayed")
            return

#    def SampleFinalTauParams(self):
        #Should add more parameters here to alter energy loss models
        #define parameters to pass to MMC
#        e = self.energy/units.GeV                   #MMC takes initial energy in GeV
#        den = self.GetCurrentDensity()*(units.cm**3)/units.gr   #convert density back to normal (not natural) units 
#        mult = den/2.7                              #This factor is passed to MMC which alters its default medium density paramater
        #this is the command that MMC gets
	    # -lpm turns on the LPM suppression of brems at high E. 
	    # -bs is the bremstrahlung model 
	    # -ph is the photonuclear model (Soyez)
	    # -bb 
#        cmd = 'awk "BEGIN {{for(i=0; i<1; i++) print {ein}, 100000000 }}" | time {Dir}../MMC/MMC/ammc -run -frejus -tau -medi="Frejus rock" -radius=1e6 -vcut=1.e-3 -rho={rho} -scat -lpm -bs=1 -ph=3 -bb=2 -sh=1 -ebig=1e16 -seed=1223 -tdir={Dir}../MMC/tables/'.format(ein=e, rho=mult, Dir=base_path)
#        #run MMC
#        out = subprocess.check_output(cmd, shell = True)
#	#parse output
#        #index = out.find('\n')
#        #e_final = float(out[:index])/1e3            #MMC returns energy in MeV
#        #distance = abs(float(out[index+2:-1]))/1e3  #MMC returns distance in meters. converting to km here.
#        e_final, distance = self.energy/1e9, 0
#	print('input energy')
#	print(np.log10(e))
#	print(mult)
#        cmd = 'awk "BEGIN {{for(i=0; i<1; i++) print {ein}, 1000000000 }}" | {Dir}../MMC/MMC/ammc -run -frejus -tau -medi="Frejus rock" -radius=1e6 -vcut=1.e-3 -rho={rho} -scat -lpm -bs=1 -ph=3 -bb=2 -sh=1 -ebig=1e16 -seed=1223 -tdir={Dir}../MMC/tables/'.format(ein=e, rho=mult, Dir=base_path)
        #run MMC
#        out = os.popen(cmd).read()
        #print(out)
	#parse output
#        index = out.find('\n')
#        e_final = float(out[:index])/1e3            #MMC returns energy in MeV
#        print('tau mmc')
#	print(np.log10(e_final))
#	distance = abs(float(out[index+2:-1]))/1e3  #MMC returns distance in meters. converting to km here.


#        return(e_final, distance)

    def InteractParticle(self, interaction):
        if self.particle_id == "tau_neutrino":
            
	    dNdEle = lambda y: DifferentialOutGoingLeptonDistribution(self.energy/units.GeV,self.energy*y/units.GeV,interaction = interaction)
            NeutrinoInteractionWeights = map(dNdEle,yy)
            NeutrinoInteractionWeights = NeutrinoInteractionWeights/np.sum(NeutrinoInteractionWeights)
            self.energy = self.energy*rand.choice(yy, p=NeutrinoInteractionWeights)

            if interaction == nsq.NeutrinoCrossSections_Current.CC:
                self.particle_id = "tau"
                self.SetParticleProperties()
                #self.history.append("Neutrino Interacted via CC")
                self.nCC += 1

            elif interaction == nsq.NeutrinoCrossSections_Current.NC:
                self.particle_id = "tau_neutrino"
                self.SetParticleProperties()
                self.nNC += 1
                #self.history.append("Neutrino Interacted via NC")
            return

        elif self.particle_id == "tau":
            Efin, Ladv = self.SampleFinalTauParams()
            Ladv *= units.km
            self.energy = Efin*units.GeV
            #Handle boundary cases, but don't scale by density anywhere
            has_exited = False
            while (self.region_distance + Ladv >= self.GetRegionDistance()) and (has_exited==False):
                    try:
                        Ladv -= (self.GetRegionDistance() - self.region_distance)
                        self.position += self.GetRegionDistance() - self.region_distance
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


def RollDice(initial_neutrino_energy,
             incoming_angle, index):
    region_distances, regions = GetDistancesPerSection(incoming_angle, earth_model_radii)
    densities = []
    for i in range(len(region_distances)):
        region_distances[i] = region_distances[i] * units.km
        densities.append(earth_model_densities[regions[i]] * units.gr/(units.cm**3))
    TotalDistance = np.sum(region_distances)

    FirstEvent = CasinoEvent("tau_neutrino",initial_neutrino_energy,0.0,region_distances,densities, index)
    EventCollection = [FirstEvent]

    CCcounter, NCcounter = 0,0
    while(not np.any(map(lambda e: (e.position >= TotalDistance) or (e.energy <= e.GetMass()) or (e.isCC), EventCollection))):
        for event in EventCollection:
	    if (event.isCC):
		print("WHY THE FUCK AM I HERE")
            if(event.particle_id == 'tau_neutrino'):
                #Determine how far you're going
                p1 = rand.random_sample()
                DistanceStep = event.GetProposedDistanceStep(event.GetCurrentDensity(), p1)
                #Handle boundary conditions and rescale according to column density
                has_exited = False
                while (event.region_distance + DistanceStep >= event.GetRegionDistance()) and (has_exited==False):
                    try:
                        DistanceStep -= (event.GetRegionDistance() - event.region_distance)
                        DistanceStep *= event.GetDensityRatio()
                        event.position += event.GetRegionDistance() - event.region_distance
                        event.region_distance = 0.
                        event.region_index += 1
                    except:
                        has_exited = True
                if has_exited:
                    event.position = TotalDistance
                    continue
                density = event.GetCurrentDensity()

                #now pick an interaction
                p2 = rand.random_sample()
                CC_lint = event.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.CC)
                p_int_CC = event.GetTotalInteractionLength(density) / CC_lint
                if(p2 <= p_int_CC):
                    event.position += DistanceStep
                    event.region_distance += DistanceStep
                #    event.position = TotalDistance
		    event.isCC = True
                    event.particle_id = "tau"
                    event.SetParticleProperties()
		 #   CCcounter += 1
		    #event.InteractParticle(nsq.NeutrinoCrossSections_Current.CC)
#		    print('charged current')
                else:
                    event.InteractParticle(nsq.NeutrinoCrossSections_Current.NC)
                    event.position += DistanceStep
                    event.region_distance += DistanceStep
		 #   NCcounter += 1
#		    print('neutral current')
#		    print(np.log10(event.energy/units.GeV))
		if(event.isCC):
		    continue
	    elif(event.particle_id == 'tau'):
                event.InteractParticle(interaction = nsq.NeutrinoCrossSections_Current.NC)
#		print('tau propagated')
#		print(np.log10(event.energy/units.GeV))
    #print("NC Counter: {}".format(NCcounter))
    #print("CC counter: {}".format(CCcounter))
    return EventCollection

cc_left = True
propagated_stack = []
inds_left = range(nevents)

taus_e = []
nus_e = []
counter = 0
iter_energies = eini
while cc_left:
    print(counter)
    counter+=1
    cc_stack = []
    for j in range(len(inds_left) - 1, -1, -1):
	#print(j)
        i = inds_left[j]
	#print(i)
        out = RollDice(iter_energies[i], thetas[i], i)[0]
	#print(out.isCC)        
        if (out.isCC):
            cc_stack.append(out)
        else:
            ind = out.index
            if (out.particle_id == 'tau'):
                taus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
                del inds_left[j]
            else:
                nus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
                del inds_left[j]
    if (len(cc_stack) > 0):
	print(len(cc_stack))
        #prop_params = [obj.energy for obj in cc_stack]
        EventCollection = DoAllCCThings(cc_stack)
        for event in EventCollection:
	    if (event.position >= event.TotalDistance):
                taus_e.append((eini[event.index], event.energy, thetas[event.index], cdf_indices[event.index]))
		delinds = np.argwhere(np.asarray(inds_left) == event.index)[0]
		print(delinds)
		print(event.index)
		print(inds_left)
                del inds_left[delinds[0]]
	    else:
                event.DecayParticle()
		iter_energies[event.index] = event.energy
    else: 
        cc_left = False


'''
#print('initial energy')
#print(np.log10(eini/units.GeV))
CasinoGame = np.array([RollDice(e, theta)[0] for (e, theta) in zip(eini, thetas)])
taus_e = []
nus_e = []

if(isgzk):
    for i, event in enumerate(CasinoGame):
        if (event.particle_id == 'tau'):
	    taus_e.append((eini[i], event.energy, thetas[i], cdf_indices[i]))
	else:
	    nus_e.append((eini[i], event.energy, thetas[i], cdf_indices[i]))

else:
    for i, event in enumerate(CasinoGame):
        if (event.particle_id == 'tau'):
	    taus_e.append((eini[i], event.energy, thetas[i]))
	else:
	    nus_e.append((eini[i], event.energy, thetas[i]))

'''
print('taus')
print(taus_e)
print('nus')
print(np.asarray(nus_e)/1e9)
#np.save('/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/taus/small_stats_taus_cosmogenic_'+str(seed)+'.npy', taus_e)
#np.save('/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/nus/small_stats_nus_cosmogenic_'+str(seed)+'.npy', nus_e)
