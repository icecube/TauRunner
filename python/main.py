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
parser.add_argument('-s',dest='seed',type=int,
    help='just an integer seed to help with output file names')
parser.add_argument('-n', dest='nevents', type=float, 
    help='how many events do you want?')
parser.add_argument('-e', dest='energy', type=float, 
    help='if you want to simulate a specific energy, pass it here in GeV')
parser.add_argument('-t', dest='theta', type=float, 
    help='nadir angle in degrees (0 is through the core)')
parser.add_argument('-gzk', dest='gzk', default=None, 
    help='Pass the file containing the CDF spline to propagate a flux model')
parser.add_argument('-spectrum', dest='spectrum', default=None, type=float, 
    help='If you want a power law, provide spectral index here')
parser.add_argument('--range', dest='range', nargs='+', 
    help='Range for injected spectrum in format "low high"')
parser.add_argument('-buff', dest='buff', type=float, default=0., 
    help="Simulate to a finite distance beneath the edge of the Earth (in km)")
parser.add_argument('-p', dest='path' , type=str, default = './', 
    help='Path to script. Default assumes you are running from within the same directory')
parser.add_argument('-d', dest='debug', default=False, action='store_true', 
    help='Do you want to print out debug statments? If so, raise this flag') 
parser.add_argument('-save', dest='save', type=str, default=None, 
    help="If saving output, provide a path here")
parser.add_argument('-onlytau', dest='onlytau', default=False, action='store_true',
    help="If you only want to save the taus, not neutrinos, raise this flag")
parser.add_argument('-water', dest='water_layer', type=float, default=0,
    help="If you you would like to add a water layer to the Earth model, enter it here in km.")
parser.add_argument('-xs', dest='xs_model', type=str, default='dipole',
    help="Enter 'CSMS' if you would like to run the simulation with a pQCD xs model")
parser.add_argument('-tau_losses', dest='tau_losses', default=True, action='store_false',
    help="Raise this flag if you want to turn off tau losses. In this case, taus will decay at rest.")

args = parser.parse_args()

save = False
isgzk = False

if ((args.seed == None) or (args.nevents == None)):
    raise RuntimeError('You must specify a seed (-s) and number of events to simulate (-n)') 
if (args.gzk == None and (args.theta==None or args.energy==None) and (args.spectrum==None)):
    raise RuntimeError('You must either pick an energy and theta, use a spectrum, or use the GZK flux')

base_path = os.path.join(args.path,'')
sys.path.append(base_path)
nevents = int(args.nevents)
seed = args.seed
debug = args.debug
xs = args.xs_model
water_layer = args.water_layer
tau_losses = args.tau_losses
if debug:
    message = ''
if args.gzk is not None:
    isgzk = True
    if not os.path.isfile(args.gzk):
        raise RuntimeError("GZK CDF Spline file does not exist")
    else:
        gzk = args.gzk
if args.save is not None:
    savedir = os.path.join(args.save, '')
    save = True
    if not os.path.isdir(savedir):
        raise RuntimeError("Directory to save output is not a valid directory")
else:
    savedir = None

from Casino import *
import Casino

cross_section_path = base_path+'../cross_sections/'

def rndm(a, b, g, size=1):
    #Random spectrum function. g is gamma+1 (use -1 for E^-2)
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

rand = np.random.RandomState(seed=seed)

if args.gzk is not None:
  # sample initial energies and incoming angles from GZK parameterization
  cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
  cdf_indices= rand.uniform(low=0., high=1.,size=nevents)
  thetas = np.arccos(cos_thetas)
  gzk_cdf = np.load(gzk, allow_pickle=True).item()
  eini = gzk_cdf(cdf_indices)*units.GeV
  if debug:
    message+="Sampled {} events from the GZK flux\n".format(nevents)
elif args.spectrum is not None:
  cdf_indices = np.ones(nevents)
  cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
  thetas = np.arccos(cos_thetas)
  eini = rndm(float(args.range[0]), float(args.range[1]), args.spectrum + 1, size=nevents)*units.GeV
  if debug:
    message+="Sampled {} events from power law\n".format(nevents)
else:
    # Use a monochromatic flux
    cdf_indices = np.ones(nevents)
    eini = np.ones(nevents)*args.energy*units.GeV
    if args.theta >= 90:
        raise ValueError("Exit angle cannot be greater than 90.")
    thetas = np.ones(nevents)*np.radians(args.theta)
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
iter_particleID = ['tau_neutrino']*nevents
iter_TauPosition = list(np.zeros(nevents))
iter_nCC = list(np.zeros(nevents))
iter_nNC = list(np.zeros(nevents))

# Run the algorithm
# All neutrinos are propagated until either exiting or undergoing a CC interaction.
# All CC interactions are handled together, and then the next iteration occurs
# This repeats until all leptons have reached the total distance
t0 = time.time()

while inds_left:
    counter += 1
    if debug:
        message+="Beginning Loop Number {}\n".format(counter)
    cc_stack = []

    for j in range(len(inds_left) - 1, -1, -1):
        i = inds_left[j] #Unique event index
    
        EventObject = CasinoEvent(iter_particleID[i],iter_energies[i], thetas[i],
			iter_positions[i], i, np.random.randint(low=1e9), iter_TauPosition[i],
			water_layer, xs_model=xs, buff = args.buff)

        out = RollDice(EventObject)

        iter_nCC[i]+=out.nCC
        iter_nNC[i]+=out.nNC      
        
        if (out.isCC):
            cc_stack.append((float(out.energy), float(out.position), int(out.index),
		 str(out.particle_id), 0, float(out.TotalDistance), float(out.GetCurrentDensity())))
            del out
        else:
            ind = int(out.index)
            if ind != i:
                message += "Index mismatch: {} {}".format(ind, i)
            if (out.particle_id == 'tau'):
                taus_e.append((eini[ind], float(out.energy), thetas[ind], cdf_indices[ind], iter_nCC[ind], iter_nNC[ind]))
               	iter_positions[out.index] = float(out.position)
                del inds_left[j]
                del out
            else:
                nus_e.append((eini[ind], float(out.energy), thetas[ind], cdf_indices[ind], iter_nCC[ind], iter_nNC[ind]))
               	iter_positions[int(out.index)] = float(out.position)
                del inds_left[j]
                del out
    if (len(cc_stack) > 0):
        if debug:
            message += "{} events passed to MMC in loop iteration {}\n".format(len(cc_stack), counter)
        EventCollection = DoAllCCThings(cc_stack, xs, tau_losses)
        for event in EventCollection:
            iter_positions[int(event[2])] = float(event[1])
            iter_energies[int(event[2])] = float(event[0])
            iter_particleID[int(event[2])] = 'tau'
            iter_TauPosition[int(event[2])] = float(event[4])
            del event

print("Simulating {} events at {} degrees took {} seconds.".format(nevents, args.theta, time.time() - t0))
nus_e = np.array(nus_e, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float), ('nCC', int), ('nNC', int)])
taus_e = np.array(taus_e, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float), ('nCC', int), ('nNC', int)])
nus_e['Theta'] *= 180. / np.pi #Give theta in degrees to user
taus_e['Theta'] *= 180. / np.pi #Give theta in degrees to user

#p_exit = 100.*(float(len(taus_e))/nevents)
#print("{} taus exited".format(p_exit))

if save:
    if isgzk:
        fluxtype = "cosmogenic"
    elif args.spectrum is not None:
        fluxtype = "powerlaw"
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
        np.save(savedir + 'nus/' + 'nus_{}_seed_{}_xs_{}_water_{}.npy'.format(fluxtype, seed, xs, water_layer), nus_e)
    np.save(savedir + 'taus/' + 'taus_{}_seed_{}_xs_{}_water_{}.npy'.format(fluxtype, seed, xs, water_layer), taus_e)
    if debug:
        print(message)
else:
    if debug:
        print(message)
    try:
        from tabulate import tabulate
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
