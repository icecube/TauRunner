#!/usr/bin/env python

import os, sys
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import argparse
import nuSQUIDSpy as nsq

info = sys.version_info
pyv  = int(info.major)

units = nsq.Const()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
tds = nsq.TauDecaySpectra()

def initialize_parser():
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
    parser.add_argument('-body', dest='body', type=str, default='earth',
        help="Raise this flag if you want to turn off tau losses. In this case, taus will decay at rest.")
    
    args = parser.parse_args()
    return args

args = initialize_parser()
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
body = args.body
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
inds_left = list(range(nevents))

taus_e = []
nus_tau_e = []
nus_electron_e = []
nus_muon_e = []
basket = []
history = []
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
			[], [], water_layer, xs_model=xs, buff=args.buff, body=body)

        out = RollDice(EventObject)

        iter_nCC[i]+=out.nCC
        iter_nNC[i]+=out.nNC   
        
        basket += out.basket
        if out.history != []: history.append((out.history, i)) # For debugging

        # Make secondaries part of the casino   
        if out.basket != []:
            SecondaryEvent = CasinoEvent(basket[len(basket)-1]["id"],basket[len(basket)-1]["energy"],thetas[i],basket[len(basket)-1]["position"],
                             i,np.random.randint(low=1e9), iter_TauPosition[i],
			                 [], [], water_layer, xs_model=xs, buff = args.buff, body=body)
            secOut = RollDice(SecondaryEvent)
        
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
                nus_tau_e.append((eini[ind], float(out.energy), thetas[ind], cdf_indices[ind], iter_nCC[ind], iter_nNC[ind]))
                if out.basket != [] and secOut is not None:
                    ind = len(basket)-1
                    if basket[ind]["id"] == "muon_neutrino":
                        nus_muon_e.append((basket[ind]["energy"], float(secOut.energy), float(secOut.position)))
                    else:
                        nus_electron_e.append((basket[ind]["energy"], float(secOut.energy), float(secOut.position)))
                    del secOut
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
nus_tau_e = np.array(nus_tau_e, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float), ('nCC', int), ('nNC', int)])
nus_muon_e = np.array(nus_muon_e, dtype = [('Eini', float), ('Eout',float), ('Finpos', float)]) # Need to tabulate
nus_electron_e = np.array(nus_electron_e, dtype = [('Eini', float), ('Eout',float), ('Finpos', float)]) # Need to tabulate
taus_e = np.array(taus_e, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float), ('nCC', int), ('nNC', int)])
basket = np.array(basket)
# history = np.array(history)
nus_tau_e['Theta'] *= 180. / np.pi #Give theta in degrees to user
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
        np.save(savedir + 'nus/' + 'nu_e_nus_{}_seed_{}_xs_{}_water_{}.npy'.format(fluxtype, seed, xs, water_layer), nus_electron_e)
        np.save(savedir + 'nus/' + 'nu_mu_nus_{}_seed_{}_xs_{}_water_{}.npy'.format(fluxtype, seed, xs, water_layer), nus_muon_e)
        np.save(savedir + 'nus/' + 'nu_tau_nus_{}_seed_{}_xs_{}_water_{}.npy'.format(fluxtype, seed, xs, water_layer), nus_tau_e)
    np.save(savedir + 'taus/' + 'taus_{}_body_{}_seed_{}_xs_{}_water_{}.npy'.format(fluxtype, body, seed, xs, water_layer), taus_e)
    if debug:
        print(message)
else:
    if debug:
        print(message)
    try:
        from tabulate import tabulate
        headers = list(nus_tau_e.dtype.names)
        nus_tau_table = tabulate(nus_tau_e, headers, tablefmt="fancy_grid")
        taus_table = tabulate(taus_e, headers, tablefmt="fancy_grid")
        print(nus_tau_table)
        print(nus_muon_e)
        print(nus_electron_e)
        print(taus_table)
        # print(basket)
        nueNumPerEvent, numuNumPerEvent = list(np.zeros(nevents)), list(np.zeros(nevents))
        for flavor, ev_Index in history:
            if flavor[0] == 'nue':
                nueNumPerEvent[ev_Index] += 1
            else: numuNumPerEvent[ev_Index] += 1
        # print(history)
        
        # print(nueNumPerEvent)
        # print(numuNumPerEvent)
        # print("")

        binsnue = []
        binsnumu = []
        for i in range(int(max(nueNumPerEvent)+1)): binsnue.append(i)
        for i in range(int(max(numuNumPerEvent)+1)): binsnumu.append(i)

        # for i in range(len(nueNumPerEvent)): nueNumPerEvent[i] = (nueNumPerEvent[i] / nevents) * 100
        # for i in range(len(numuNumPerEvent)): numuNumPerEvent[i] = (numuNumPerEvent[i] / nevents) * 100

        # print(nueNumPerEvent)
        # print(numuNumPerEvent)
        from matplotlib import pyplot as plt
        plt.hist(nueNumPerEvent, bins=binsnue, histtype='step', label=r'$\nu_e$')
        plt.hist(numuNumPerEvent, bins=binsnumu, histtype='step', label=r'$\nu_\mu$')
        # plt.yscale('log', nonposy='clip')
        plt.xlabel('# of secondaries per nu_tau')
        plt.ylabel('%')
        plt.legend(loc='best')
        plt.show()
    except ImportError:
        print("Outgoing Neutrinos: ")
        print(nus_tau_e)
        print(nus_muon_e)
        print(nus_electron_e)
        print("Outgoing Taus: ")
        print(taus_e)
        print(basket)
        # print(history)
