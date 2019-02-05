import os 
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import scipy as sp
import nuSQUIDSpy as nsq
import pickle
import argparse
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser()
parser.add_argument('-s',dest='seed',type=int,help='just an integer seed to help with output file names')
parser.add_argument('-e',dest='eini',type=float,help='initial nutau energy in GeV')
parser.add_argument('-d', dest='total_distance', type=float, help='total distance you would like to propagate the neutrino/tau')
parser.add_argument('-n', dest='nevents', type=float, help='how many events do you want?')
args = parser.parse_args()

if (not (args.seed or args.eini or args.total_distance or args.nevents)):
    raise RuntimeError('You must specify a seed (-s), an initial energy (-e), total propagation distance (-d) and number of events to simulate (-n)') 

seed = args.seed
eini = args.eini
total_distance = args.total_distance
nevents = int(args.nevents)

units = nsq.Const()
gr = nsq.GlashowResonanceCrossSection()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
tds = nsq.TauDecaySpectra()

File = open('taus_rock_ALLM97_lpm_vcut1e-3_MMC_Ordered_weighted.pickle','r')
tau_loss_array = pickle.load(File)

#cross sections from 1e8 - 1e16 GeV patched with nuSQUIDS 
#temporary until EHE cross sections are added to nuSQUIDS

######################################
charged = np.load('nuXS_CC_8-16.npy')
neutral = np.load('nuXS_NC_8-16.npy')
energies = np.logspace(17, 25, 500)
energies = (energies[:-1] + energies[1:])/2

log_e = np.log10(energies)
log_XS_CC = np.log10(charged)
log_XS_NC = np.log10(neutral)

f_NC = interp1d(log_e, log_XS_NC)
f_CC = interp1d(log_e, log_XS_CC)

#####################################

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def SampleFinalTauParams(e_in):
    #below a PeV decay on-the-spot
    if (e_in < 1e6):
        return(e_in, 0.)
    else:
        eins = np.logspace(6,16,500)
        e_in = find_nearest(eins, e_in)
        choices = tau_loss_array[e_in]['choices']
        weights = tau_loss_array[e_in]['weights']
        choice_index = np.arange(0, len(choices))

        index = np.random.choice(choice_index, p=weights)
        (e_final, distance) = choices[index]

        return(10**e_final, 10**distance)

def TotalNeutrinoCrossSection(enu,
                              flavor = nsq.NeutrinoCrossSections_NeutrinoFlavor.tau,
                              neutype = nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
                              interaction = nsq.NeutrinoCrossSections_Current.NC):

    if(np.log10(enu) > log_e[0]):

        if(interaction == nsq.NeutrinoCrossSections_Current.NC):
            return((10**f_NC(np.log10(enu)))*(units.cm)**2)
        else:
            return((10**f_CC(np.log10(enu)))*(units.cm)**2)
    else:

        return dis.TotalCrossSection(enu,flavor,neutype,interaction)*(units.cm)**2

def DifferentialOutGoingLeptonDistribution(enu_in,enu_out,
                                       flavor = nsq.NeutrinoCrossSections_NeutrinoFlavor.tau,
                                       neutype = nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
                                       interaction = nsq.NeutrinoCrossSections_Current.NC
                                    ):
    diff = dis.SingleDifferentialCrossSection(enu_in,enu_out,flavor,neutype,interaction)
    return diff

Etau = 100.
zz = np.linspace(0.0,1.0,500)[1:-1]
dNTaudz = lambda z: tds.TauDecayToAll(Etau, Etau*z)
TauDecayWeights = np.array(map(dNTaudz,zz))
TauDecayWeights = TauDecayWeights/np.sum(TauDecayWeights)

yy = np.linspace(0.0,1.0,300)[1:-1]

proton_mass = 0.938*units.GeV


class CasinoEvent(object):
    def __init__(self, particle_id, energy, position):
        ## need to add tree in an efficient way

        self.particle_id = particle_id
        self.energy = energy
        self.position = position
        self.SetParticleProperties()
        self.history = ["Event created as " + particle_id]

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

    def DecayParticle(self):
        if self.particle_id == "tau_neutrino":
            self.history.append("Neutrino decayed???")
            return
        if self.particle_id == "tau":
            self.energy = self.energy*np.random.choice(zz, p=TauDecayWeights)
            self.particle_id = "tau_neutrino"
            self.SetParticleProperties()
            self.history.append("Tau decayed")
            return

    def InteractParticle(self, interaction):
        if self.particle_id == "tau_neutrino":

	    if(self.energy/units.GeV > 1e12):
  		dNdEle = lambda y: DifferentialOutGoingLeptonDistribution(1e12*units.GeV,1e12*units.GeV*y,
                                                                      interaction = interaction)
	    else:
                dNdEle = lambda y: DifferentialOutGoingLeptonDistribution(self.energy,self.energy*y,
                                                                      interaction = interaction)
            NeutrinoInteractionWeights = map(dNdEle,yy)
            NeutrinoInteractionWeights = NeutrinoInteractionWeights/np.sum(NeutrinoInteractionWeights)
            self.energy = self.energy*np.random.choice(yy, p=NeutrinoInteractionWeights)

            if interaction == nsq.NeutrinoCrossSections_Current.CC:
                self.particle_id = "tau"
                self.SetParticleProperties()
		self.history.append("Neutrino Interacted via CC")

            elif interaction == nsq.NeutrinoCrossSections_Current.NC:
                self.particle_id = "tau_neutrino"
                self.SetParticleProperties()
		self.history.append("Neutrino Interacted via NC")
            return
        if self.particle_id == "tau":
            Efin, Ladv = SampleFinalTauParams(self.energy/units.GeV)
            self.energy = Efin*units.GeV
            self.position += Ladv*units.km
            self.history.append("Tau Interacted")
	    self.history.append(" moved "+str(Ladv)+" km and has "+str(Efin)+" GeV left") 
            self.DecayParticle()
            return

    def PrintParticleProperties(self):
        print "id", self.particle_id, \
              "energy ", self.energy/units.GeV, " GeV", \
              "position ", self.position/units.km, " km"


def RollDice(initial_neutrino_energy,
             TotalDistance,
             density = 2.6*units.gr/(units.cm**3)):

    FirstEvent = CasinoEvent("tau_neutrino",initial_neutrino_energy,0.0)
    EventCollection = [FirstEvent]

    while(not np.any(map(lambda e: (e.position >= TotalDistance) or (e.energy <= e.GetMass()), EventCollection))):
        for event in EventCollection:

            if(event.particle_id == 'tau_neutrino'):

                #Determine how far you're going
                p1 = np.random.random_sample()
                DistanceStep = event.GetProposedDistanceStep(density, p1)

                #Check to see if your particle escapes before interacting
                if(event.position + DistanceStep >= TotalDistance):
                    event.position = TotalDistance
                    continue

                #now pick an interaction
                p2 = np.random.random_sample()
                CC_lint = event.GetInteractionLength(density, interaction=nsq.NeutrinoCrossSections_Current.CC)
                p_int_CC = event.GetTotalInteractionLength(density) / CC_lint

                if(p2 <= p_int_CC):
                    event.InteractParticle(nsq.NeutrinoCrossSections_Current.CC)
                    event.position += DistanceStep
                else:
                    event.InteractParticle(nsq.NeutrinoCrossSections_Current.NC)
                    event.position += DistanceStep

	    elif(event.particle_id == 'tau'):
                event.InteractParticle(interaction = nsq.NeutrinoCrossSections_Current.NC)

    return EventCollection

eini = eini*units.GeV
print(np.log10(eini))
total_distance = total_distance*units.km
CasinoGame = np.array([RollDice(eini, total_distance)[0] for i in xrange(nevents)])

taus_e = []
nus_e = []

for event in CasinoGame:
    if (event.particle_id == 'tau'):
        taus_e.append(event.energy)
    else:
        nus_e.append(event.energy)

np.save('./taus/taus_onTheOtherSide'+str(seed)+'_'+str(eini/units.GeV)+'GeV.npy', taus_e)
np.save('./nus/nus_onTheOtherSide'+str(seed)+'_'+str(eini/units.GeV)+'GeV.npy', nus_e)
