#!/usr/bin/env python

import os, sys
import json
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import argparse

from taurunner.modules import units, make_outdir, todaystr, cleanup_outdir
from taurunner.track import Chord
from taurunner.Casino import *

def initialize_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',dest='seed',type=int,
        help='just an integer seed to help with output file names')
    parser.add_argument('-n', dest='nevents', type=float, 
        help='how many events do you want?')
    parser.add_argument('-flavor', dest='flavor', type=int, default=3,
        help='neutrino flavor (default is nutau): 1 for nue 2 for numu 3 for nutau')
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
    parser.add_argument('-d', dest='debug', default=False, action='store_true', 
        help='Do you want to print out debug statments? If so, raise this flag') 
    parser.add_argument('-save', dest='save', type=str, default=None, 
        help="If saving output, provide a path here")
    parser.add_argument('-water', dest='water_layer', type=float, default=0,
        help="If you you would like to add a water layer to the Earth model, enter it here in km.")
    parser.add_argument('-xs', dest='xs_model', type=str, default='dipole',
        help="Enter 'CSMS' if you would like to run the simulation with a pQCD xs model")
    parser.add_argument('-losses', dest='losses', default=True, action='store_false',
        help="Raise this flag if you want to turn off tau losses. In this case, taus will decay at rest.")
    parser.add_argument('--body', dest='body', type=str, default='earth',
        help="Raise this flag if you want to turn off tau losses. In this case, taus will decay at rest.")
    parser.add_argument('--depth', dest='depth', type=float, default=0.0,
        help="Depth of the detector in km.")
    
    args = parser.parse_args()
    return args

def propagate_neutrinos(nevents, seed, flavor=3, energy=None, theta=None,
    gzk=None, spectrum=None, e_range=" ", debug=False, save=None, onlytau=False,
    water_layer=0., xs_model='dipole', losses=True, body='earth', depth=0., 
    return_res=True):
    r'''
    Main simulation code. Propagates a flux of neutrinos and returns or
    saves the outgoing particles

    Parameters
    ----------
    nevents: int
        Number of neutrinos to propagate
    seed: int
        Random number seed
    flavor: int, default=3
        Neutrino flavor. select 1, 2, 3, for nue, numu and nutau respectively.
    energy: float, default=None
        Energy in GeV for injecting a monochromatic flux of neutrinos
    theta: float, default=None
        nadir angle in degrees for fixed angle simulation
    gzk: str, default=None
        Path to a gzk cdf spline to simulate a specific flux
    spectrum: float, default=None
        If not None, inject flux as a power-law with given spectral index
    e_range: str, default=" "
        If injecting power-law, the energy range to use, format is 'low_en high_en'
    debug: bool, default=False
        If true, print various messages as the simulation runs
    save: str, default=None
        If not None, a path for saving output
    onlytau: bool, default=False
        If true, only save the outgoing taus, not the neutrinos
    water_layer: float, default=0.
        Add a layer of water on top of the body
    xs_model: str, default='dipole'
        Cross section model to use
    losses: bool, default=True
        If false, do not include stochastic losses
    body: str, default='earth'
        Name of the body which we propagate through
    buff: float, default=0.
        Buffer in km between the end of propagation and the detector
    return_res: bool, default=True
        If true, return the results
    
    Returns:
        np.recarray of the outgoing leptons
    
    '''

    if nevents is None:
        raise RuntimeError('You must specify a number of events to simulate (-n)') 
    if (gzk == None and theta == None) or (energy ==None and spectrum ==None):
        raise RuntimeError('You must either pick an energy and theta, use a spectrum, or use the GZK flux')
   
    if seed is None:
        seed = int(float(savedir.split('/')[-1].replace('_', ''))) % 2**32
    else:
        seed = seed

    print('Beggining simulation')
    nevents     = int(nevents)
    depth       = depth*units.km
    gzk         = gzk
    theta       = theta

    if(body=='earth'):
        from taurunner.body import Earth
        body = Earth
    elif(body=='sun'):
        from taurunner.body import HZ_Sun
        body = HZ_Sun

    if debug:
        message = ''

    def rndm(a, b, g, size=1):
        #Random spectrum function. g is gamma+1 (use -1 for E^-2)
        r = np.random.random(size=size)
        if g == 0.:
            # E^-1 is uniform sampling in log space
            log_es = (np.log10(b) - np.log10(a)) * r + np.log10(a)
            return 10.**log_es
        ag, bg = a**g, b**g
        return (ag + (bg - ag)*r)**(1./g)

    rand = np.random.RandomState(seed=seed)

    if gzk is not None:
        if not os.path.isfile(gzk):
            raise RuntimeError("GZK CDF Spline file does not exist")
        # sample initial energies and incoming angles from GZK parameterization
        cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
        cdf_indices= rand.uniform(low=0., high=1.,size=nevents)
        thetas = np.arccos(cos_thetas)
        gzk_cdf = np.load(gzk, allow_pickle=True).item()
        eini = gzk_cdf(cdf_indices)*units.GeV
        if debug:
            message+="Sampled {} events from the GZK flux\n".format(nevents)
    elif spectrum is not None:
        cdf_indices = np.ones(nevents)
        if theta is None:
            cos_thetas = rand.uniform(low=0., high=1.,size=nevents)
            thetas = np.arccos(cos_thetas)
        else:
            if theta >= 90:
                raise ValueError("Exit angle cannot be greater than 90.")
            if theta is not None:
                thetas = np.ones(nevents)*np.radians(theta)
            else:
                cos_thetas = rand.uniform(low=0., high=1., size=nevents)
                thetas = np.arccos(cos_thetas)
        eini = rndm(float(e_range[0]), float(e_range[1]), spectrum + 1, size=nevents)*units.GeV
        if debug:
            message+="Sampled {} events from power law\n".format(nevents)
    else:
        # Use a monochromatic flux
        cdf_indices = np.ones(nevents)
        eini = np.ones(nevents)*energy*units.GeV
        if theta is not None:
            if theta >= 90:
                raise ValueError("Exit angle cannot be greater than 90.")
            thetas = np.ones(nevents)*np.radians(theta)
        else:
            cos_thetas = rand.uniform(low=0., high=1., size=nevents)
            thetas = np.arccos(cos_thetas)
        if debug:
            message+="Sampled {} events from monochromatic flux\n".format(nevents)

    cc_left = True
    propagated_stack = []
    inds_left = list(range(nevents))

    output  = []
    counter = 0
    iter_energies = list(eini)[:]
    iter_positions = list(np.zeros(nevents))
    if(flavor==2):
        iter_particleID = np.ones(nevents, dtype=int)*14
        flavors = ['mu']*nevents
    elif(flavor==3):
        iter_particleID = np.ones(nevents, dtype=int)*16
        flavors = ['tau']*nevents
    iter_ChargedPosition = list(np.zeros(nevents))
    iter_nCC = list(np.zeros(nevents))
    iter_nNC = list(np.zeros(nevents))

    # Run the algorithm
    # All neutrinos are propagated until either exiting or undergoing a CC interaction.
    # All CC interactions are handled together, and then the next iteration occurs
    # This repeats until all leptons have reached the total distance
    t0 = time.time()
    tracks  = {theta:Chord(theta=theta, depth=depth/body.radius) for theta in set(thetas)}
    while inds_left:
        counter += 1
        if debug:
            message+="Beginning Loop Number {}\n".format(counter)
        cc_stack = []

        for j in range(len(inds_left) - 1, -1, -1):
            i = inds_left[j] #Unique event index

            particle = Particle(iter_particleID[i], flavors[i], iter_energies[i], 
                                thetas[i], iter_positions[i], i, rand.randint(low=1e9),
                                iter_ChargedPosition[i], xs_model=xs_model)
            my_track = tracks[thetas[i]]
            out = Propagate(particle, my_track, body)

            iter_nCC[i]+=out.nCC
            iter_nNC[i]+=out.nNC      
            if (out.survived==False):
                # these were absorbed. we record them in the output with outgoing energy 0
                output.append((eini[ind], 0., thetas[ind], cdf_indices[ind], iter_nCC[ind], iter_nNC[ind], out.ID))
                iter_positions[int(out.index)] = float(out.position)
                del inds_left[j]
                del out
            elif (out.isCC):
                current_distance=my_track.x_to_d(out.position)*body.radius
                current_x = out.position
                total_distance=my_track.x_to_d(1.)*body.radius
                current_density=body.get_density(my_track.x_to_r(out.position))
                cc_stack.append((float(out.energy), current_x, float(current_distance), int(out.index),
                 int(out.ID), 0, float(total_distance), float(current_density), str(out.flavor)))
                del out
            else:
                ind = int(out.index)
                if ind != i:
                    message += "Index mismatch: {} {}".format(ind, i)
                    raise RuntimeError('Index mismatch -- particles are getting jumbled somewhere (thats a bad thing)')
                output.append((eini[ind], float(out.energy), thetas[ind], cdf_indices[ind], iter_nCC[ind], iter_nNC[ind], out.ID))
                iter_positions[int(out.index)] = float(out.position)
                del inds_left[j]
                del out
        if (len(cc_stack) > 0):
            if debug:
                message += "{} events passed to MMC in loop iteration {}\n".format(len(cc_stack), counter)
            EventCollection = DoAllCCThings(cc_stack, xs_model, flavor, losses)
            for event in EventCollection:
                iter_positions[int(event[3])] = float(event[1])
                iter_energies[int(event[3])] = float(event[0])
                iter_particleID[int(event[3])] = int(event[4])
                iter_ChargedPosition[int(event[3])] = float(event[5])
                del event

    print("Simulating {} events at {} degrees took {} seconds.".format(nevents, theta, time.time() - t0))

    output = np.array(output, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float), ('nCC', int), ('nNC', int), ('PDG_Encoding', int)])
    output['Theta'] *= 180. / np.pi #Give theta in degrees to user

    return output

if __name__ == "__main__":

    args = initialize_parser()

    if args.save is not None:
        savedir = os.path.join(args.save, '')
        save    = True
        if not os.path.isdir(savedir):
            raise RuntimeError("Directory to save output is not a valid directory")
        savedir = make_outdir(savedir, todaystr)
        os.mkdir(savedir)
        params_file = savedir+"/params.json"
        output_file = savedir+'/output.npy'

        d = vars(args)
        # Check this
        d['seed'] = args.seed
        j = json.dumps(d)
        f = open(params_file,"w")
        f.write(j)
        f.close()

    try:
        result = propagate_neutrinos(args.nevents, args.seed, flavor=args.flavor, 
            energy=args.energy, theta=args.theta, gzk=args.gzk, 
            spectrum=args.spectrum, e_range=args.range, debug=args.debug,
            save=args.save, water_layer=args.water_layer, xs_model=args.xs_model,
            losses=args.losses, body=args.body, depth=args.depth, return_res=False)

        if args.save:
            if args.gzk is not None:
                fluxtype = "cosmogenic"
            elif args.spectrum is not None:
                fluxtype = "powerlaw"
            else:
                fluxtype = "monochromatic_{}_{}".format(args.energy, args.theta)
            np.save(output_file, result)
            if args.debug:
                print(message)
        else:
            if args.debug:
                print(message)
            try:
                from tabulate import tabulate
                headers = list(result.dtype.names)
                out_table = tabulate(result, headers, tablefmt="fancy_grid")
                print(out_table)
            except ImportError:
                print("Outgoing Particles: ")
                print(result)
    
     except BaseException as err:
         cleanup_outdir(savedir, output_file, params_file)
         raise err
 
