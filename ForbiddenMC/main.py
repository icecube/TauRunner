#!/usr/bin/env python

import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import os, sys
import argparse
import nuSQUIDSpy as nsq

units = nsq.Const()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
tds = nsq.TauDecaySpectra()

parser = argparse.ArgumentParser()
parser.add_argument('-s',dest='seed',type=int,help='just an integer seed to help with output file names')
parser.add_argument('-n', dest='nevents', type=float, help='how many events do you want?')
parser.add_argument('-gzk', dest='gzk', default=False, action='store_true', help='do you want to propagate the GZK flux? if so, raise this flag, and raise your flag, and raise your flag, and raise it.. when i get older...')
parser.add_argument('-e', dest='energy', type=float, help='if you want to simulate a specific energy, pass it here in GeV')
parser.add_argument('-t', dest='theta', type=float, help='zenith angle in radians where 0 is through the core')
parser.add_argument('-p', dest='path' , type=str, help='path to script')
parser.add_argument('-d', dest='debug', default=False, action='store_true', help='Do you want to print out debug statments? If so, raise this flag') 
parser.add_argument('-save', dest='save', default=False, action='store_true', help='Do you want to save the output or print the output? If save, raise this flag')
parser.add_argument('-savedir', dest='savedir', type=str, default='./', help="If saving output, where would you like to save output?")
parser.add_argument('-onlytau', dest='onlytau', default=False, action='store_true', help='Only save the taus and not neutrinos if this flag is raised')
parser.add_argument('spectrum', dest='spectrum', default=False, action='store_true', help='Inject an E^-2 spectrum in a specific bin')
parser.add_argument('spectral_index', dest='spectral_index', default=-2., type=float, help='Spectral index within the bin')
parser.add_argument('low', dest='low', type=float, default=1e3, help='lower bound if using an injected spectrum')
parser.add_argument('high', dest='high', type=float, default=1e6, help='upper bound if using an injected spectrum')
args = parser.parse_args()

if ((args.seed == None) or (args.nevents == None)):
    raise RuntimeError('You must specify a seed (-s) and number of events to simulate (-n)') 
if (not (args.gzk) and (args.theta==None and args.energy==None)):
    raise RuntimeError('You must either pick an energy and theta or use the GZK flux, bud')

base_path = os.path.join(args.path,'')
sys.path.append(base_path)
#print(sys.path[-1])

from Casino import *
import Casino
#from CrossSections import *
#import CrossSections

seed = args.seed
debug = args.debug
if debug:
    message = ''
nevents = int(args.nevents)
isgzk = args.gzk
save = args.save
savedir = os.path.join(args.savedir, '')

#base_path = '/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/'
cross_section_path = base_path+'../cross_sections/'

#Initialize cross section class
#CrossSection = CrossSections.CrossSections(cross_section_path, seed)
#if debug:
#    message+="Initialized Cross Secions\n"

def rndm(a, b, g, size=1):
    #Random spectrum function. g is gamma+1 (use -1 for E^-2)
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

rand = np.random.RandomState(seed=seed)

if(isgzk):
  # sample initial energies and incoming angles from GZK parameterization
  cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
  thetas = np.arccos(cos_thetas)
  #thetas = np.zeros(nevents)

  gzk_cdf = np.load(base_path+'gzk_cdf_phi_spline.npy').item()
  cdf_indices = rand.uniform(size=nevents)
  eini = gzk_cdf(cdf_indices)*units.GeV
  #eini = np.ones(nevents)*1e9*units.GeV
  if debug:
    message+="Sampled {} events from the GZK flux\n".format(nevents)
elif spectrum:
  cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
  thetas = np.arccos(cos_thetas)
  eini = rndm(args.low, args.high, args.spectral_index + 1, size=nevents)
else:
    # Use a monochromatic flux
    cdf_indices = np.ones(nevents)
    eini = np.ones(nevents)*args.energy*units.GeV
    thetas = np.ones(nevents)*args.theta
    if debug:
        message+="Sampled {} events from monochromatic flux\n".format(nevents)

cc_left = True
propagated_stack = []
inds_left = range(nevents)

taus_e = []
nus_e = []
counter = 0
iter_energies = list(eini)[:]
iter_positions = list(np.zeros(nevents))

# Run the algorithm
while cc_left:
    counter += 1
    if debug:
        message+="Beginning Loop Number {}\n".format(counter)
    cc_stack = []
    low_en_cc = []
    for j in range(len(inds_left) - 1, -1, -1):
    #for i in inds_left[::-1]:
        i = inds_left[j]
#        print("lookin at {}".format(i))
        EventObject = CasinoEvent("tau_neutrino",iter_energies[i], thetas[i], iter_positions[i], i, np.random.randint(low=1e9))
        out = RollDice(EventObject)       
        if (out.isCC):
            if(out.energy/units.GeV <= 1e5):
                out.DecayParticle()
                iter_positions[out.index] = out.position
                iter_energies[out.index] = out.energy 
                low_en_cc.append(out)
            else:
                cc_stack.append(out)
        else:
            ind = out.index
            if ind != i:
                print ("THIS IS A WRONG INDEX {} {}".format(ind, i))
            if (out.particle_id == 'tau'):
                taus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
               	iter_positions[out.index] = out.position
                del inds_left[j]
                del out
            else:
                if not args.onlytau:
                    nus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
               	iter_positions[out.index] = out.position
                del inds_left[j]
                del out
    if (len(cc_stack) > 0):
        if debug:
            message += "{} events passed to MMC in loop iteration {}\n".format(len(cc_stack), counter)
            #print("{} events passed to MMC in loop iteration {}\n".format(len(cc_stack), counter))
        EventCollection = DoAllCCThings(cc_stack)
        for event in EventCollection:
            if (event.position >= event.TotalDistance):
                taus_e.append((eini[event.index], event.energy, thetas[event.index], cdf_indices[event.index]))
                delinds = np.argwhere(np.asarray(inds_left) == event.index)[0]
                iter_positions[event.index] = event.position
                del inds_left[delinds[0]]
                del event
            else:
                event.DecayParticle()
                iter_positions[event.index] = event.position
                iter_energies[event.index] = event.energy
                del event
    else: 
        if len(low_en_cc) == 0:
            cc_left = False

#print("INDICES LEFT: {}".format(inds_left))
nus_e = np.array(nus_e, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float)])
taus_e = np.array(taus_e, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float)])
#print(np.array(iter_positions) / units.km)

if save:
    if isgzk:
        fluxtype = "cosmogenic"
    else:
        fluxtype = "monochromatic_{}_{}".format(args.energy, args.theta)
    try:
        os.mkdir(savedir + 'nus/')
        os.mkdir(savedir + 'taus/')
        if debug:
            message += "Created subdirectories for nus and taus\n"
    except OSError as e:
        if debug:
            message += "Subdirectories already existed\n"
    if not args.onlytau:
        np.save(savedir + 'nus/' + 'nus_{}_seed_{}.npy'.format(fluxtype, seed), nus_e)
    np.save(savedir + 'taus/' + 'taus_{}_seed_{}.npy'.format(fluxtype, seed), taus_e)
    if debug:
        print(message)
else:
    if debug:
        print(message)
    try:
        from tabulate import tabulate
        #headers = ["Energy In", "Energy Out", "Theta", "CDF Index"]
        headers = list(nus_e.dtype.names)
        nus_table = tabulate(nus_e, headers, tablefmt="fancy_grid")
        taus_table = tabulate(taus_e, headers, tablefmt="fancy_grid")
        print(nus_table)
        print(taus_table)
    except ImportError:
        print("Outgoing Neutrinos: ")
        print(nus_e)
        print("Outgoing Taus: ")
        print(taus_e)

