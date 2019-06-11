#!usr/bin/env python

from Casino import *
import Casino
import CrossSections
import argparse, os
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
args = parser.parse_args()

if (not (args.seed or args.nevents)):
    raise RuntimeError('You must specify a seed (-s) and number of events to simulate (-n)') 
if (not (args.gzk) and not (args.theta and args.energy)):
    raise RuntimeError('You must either pick an energy and theta or use the GZK flux, bud')

seed = args.seed
debug = args.debug
if debug:
    message = ''
nevents = int(args.nevents)
isgzk = args.gzk
base_path = args.path
save = args.save
savedir = os.path.join(args.savedir, '')
#base_path = '/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/'
cross_section_path = base_path+'../cross_sections/'

#Initialize cross section class
CrossSection = CrossSections.CrossSections(cross_section_path, seed)
if debug:
    message+="Initialized Cross Secions\n"

rand = np.random.RandomState(seed=seed)

if(isgzk):
  # sample initial energies and incoming angles from GZK parameterization
  cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
  thetas = np.arccos(cos_thetas)
  #thetas = np.zeros(nevents)

  gzk_cdf = np.load(base_path+'gzk_cdf_phi_spline.npy').item()
  cdf_indices = CrossSection.rand.uniform(size=nevents)
  eini = gzk_cdf(cdf_indices)*units.GeV
  #eini = np.ones(nevents)*1e9*units.GeV
  if debug:
    message+="Sampled {} events from the GZK flux\n".format(nevents)
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
iter_energies = eini
iter_positions = list(np.zeros(nevents))

# Run the algorithm
while cc_left:
    counter += 1
    if debug:
        message+="Beginning Loop Number {}\n".format(counter)
    cc_stack = []
    low_en_cc = []
    for j in range(len(inds_left) - 1, -1, -1):
        i = inds_left[j]
        EventObject = CasinoEvent("tau_neutrino",iter_energies[i], thetas[i], iter_positions[i], i, CrossSection, seed+i)
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
            if (out.particle_id == 'tau'):
                taus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
               	iter_positions[out.index] = out.position
                del inds_left[j]
            else:
                nus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
               	iter_positions[out.index] = out.position
                del inds_left[j]
    if (len(cc_stack) > 0):
        if debug:
            message += "{} events passed to MMC in loop iteration {}\n".format(len(cc_stack), counter)
        EventCollection = DoAllCCThings(cc_stack)
        for event in EventCollection:
            if (event.position >= event.TotalDistance):
                taus_e.append((eini[event.index], event.energy, thetas[event.index], cdf_indices[event.index]))
                delinds = np.argwhere(np.asarray(inds_left) == event.index)[0]
                iter_positions[event.index] = event.position
                del inds_left[delinds[0]]
            else:
                event.DecayParticle()
                iter_positions[event.index] = event.position
                iter_energies[event.index] = event.energy
    else: 
        if len(low_en_cc) == 0:
            cc_left = False

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
    np.save(savedir + 'nus/' + 'nus_{}_seed_{}.npy'.format(fluxtype, seed), nus_e)
    np.save(savedir + 'taus/' + 'taus_{}_seed_{}.npy'.format(fluxtype, seed), taus_e)
    if debug:
        print(message)
else:
    if debug:
        print(message)
    try:
        from tabulate import tabulate
        headers = ["Energy In", "Energy Out", "Theta", "CDF Index"]
        nus_table = tabulate(nus_e, headers, tablefmt="fancy_grid")
        taus_table = tabulate(taus_e, headers, tablefmt="fancy_grid")
        print(nus_table)
        print(taus_table)
    except ImportError:
        print("Outgoing Neutrinos: ")
        print(nus_e)
        print("Outgoing Taus: ")
        print(taus_e)

