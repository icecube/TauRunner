import numpy as np
from taurunner.utils import units
from scipy.interpolate import InterpolatedUnivariateSpline as iuvs
from importlib.resources import path
from taurunner.resources import secondaries_splines

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

proton_mass = ((0.9382720813+0.9395654133)/2.)*units.GeV
NeutrinoDifferentialEnergyFractions = np.linspace(0.0,1.0,1000)[1:-1]
TauDecayFractions = np.linspace(0.0,1.0,1000)[1:-1]
Etau = 100.
dNTaudz = lambda z: TauDecayToAll(Etau, Etau*z, -1.)
TauDecayWeights = np.array(list(map(dNTaudz,TauDecayFractions)))
TauDecayWeights = np.divide(TauDecayWeights, np.sum(TauDecayWeights))

#########################################################
#### SECONDARIES ########################################
#########################################################

with path(secondaries_splines, 'secondaries_nuebar_spline.npy') as p:
        nue_path = str(p)
with path(secondaries_splines, 'secondaries_numubar_spline.npy') as p:
        numu_path = str(p)

antinue_cdf = np.load(nue_path, allow_pickle=True).item()
antinumu_cdf = np.load(numu_path, allow_pickle=True).item()

def SampleSecondariesEnergyFraction(u, cdf):
    # check if random sample u is in the range where the spline is defined
    if(u<1e-3):
        return 10**cdf(-3)
    else:
        return 10**cdf(np.log10(u))
