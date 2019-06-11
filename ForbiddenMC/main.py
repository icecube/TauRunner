#!usr/bin/env python

from Casino import *
import Casino
import CrossSections
import argparse
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
args = parser.parse_args()

if (not (args.seed or args.nevents)):
    raise RuntimeError('You must specify a seed (-s) and number of events to simulate (-n)') 

seed = args.seed
nevents = int(args.nevents)
isgzk = args.gzk
base_path = args.path
#base_path = '/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/'
cross_section_path = base_path+'../cross_sections/'

#Initialize cross section class
CrossSection = CrossSections.CrossSections(cross_section_path, seed)

if(isgzk):
  # sample initial energies and incoming angles
  cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
  thetas = np.arccos(cos_thetas)
  #thetas = np.zeros(nevents)

  gzk_cdf = np.load(base_path+'gzk_cdf_phi_spline.npy').item()
  cdf_indices = CrossSection.rand.uniform(size=nevents)
  eini = gzk_cdf(cdf_indices)*units.GeV
  #eini = np.ones(nevents)*1e9*units.GeV
else:
    eini = np.ones(nevents)*args.energy*units.GeV
    thetas = np.ones(nevents)*args.theta

cc_left = True
propagated_stack = []
inds_left = range(nevents)

taus_e = []
nus_e = []
counter = 0
iter_energies = eini
iter_positions = list(np.zeros(nevents))
#print(iter_positions)
#print(inds_left)

while cc_left:
    print('\n \n')
    counter+=1
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
		print(out.index)
		low_en_cc.append(out)
	    else:
		cc_stack.append(out)
        else:
            ind = out.index
            if (out.particle_id == 'tau'):
                taus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
               	iter_positions[out.index] = out.position
		#if(iter_positions[out.index] > out.TotalDistance):
		#    print(out.energy, out.nCC, out.index, out.nNC, out.isCC)
	        #print('deleting index: {}'.format(inds_left[j]))
		del inds_left[j]
            else:
                nus_e.append((eini[ind], out.energy, thetas[ind], cdf_indices[ind]))
               	iter_positions[out.index] = out.position
		#if(iter_positions[out.index] > out.TotalDistance):
		#    print(out.energy, out.nCC, out.index, out.nNC, out.isCC)
		#print('deleting index: {}'.format(inds_left[j]))
		del inds_left[j]
    if (len(cc_stack) > 0):
	print(len(cc_stack))
        EventCollection = DoAllCCThings(cc_stack)
        for event in EventCollection:
	    if (event.position >= event.TotalDistance):
                taus_e.append((eini[event.index], event.energy, thetas[event.index], cdf_indices[event.index]))
		delinds = np.argwhere(np.asarray(inds_left) == event.index)[0]
 #               print('deleting index: {}'.format(inds_left[delinds[0]]))
		iter_positions[event.index] = event.position
#		if(iter_positions[event.index] > event.TotalDistance):
#		    print(event.energy, event.nCC, event.index, event.nNC, event.isCC)
		del inds_left[delinds[0]]
	    else:
                event.DecayParticle()
		iter_positions[event.index] = event.position
		iter_energies[event.index] = event.energy
    else: 
        if len(low_en_cc) == 0:
	    cc_left = False

print(np.asarray(iter_positions)/units.km)

#print('initial energy')
#print(np.log10(eini/units.GeV))
#CasinoGame = np.array([RollDice(e, theta)[0] for (e, theta) in zip(eini, thetas)])
#taus_e = []
#nus_e = []
#np.save('/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/taus/small_stats_taus_cosmogenic_'+str(seed)+'.npy', taus_e)
#np.save('/data/user/isafa/ANITA/features/TauDragon/ForbiddenMC/nus/small_stats_nus_cosmogenic_'+str(seed)+'.npy', nus_e)
#
print('GOTEM')
print(nus_e)

