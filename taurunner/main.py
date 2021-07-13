#!/usr/bin/env python

import os, sys, time
import json
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import argparse

from taurunner.track import Chord
from taurunner.body import *
from taurunner.cross_sections import CrossSections
from taurunner.Casino import *
from taurunner.particle import Particle


def initialize_parser(): # pragma: no cover
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        dest='seed',
                        type=int, 
                        default=None,
                        help='just an integer seed to help with output file names'
                       )
    parser.add_argument('-n', 
                        dest='nevents', 
                        type=int, 
                        default=0,
                        help='how many events do you want?'
                       )
    parser.add_argument('-flavor', 
                        dest='flavor', 
                        type=int, default=16,
                        help='neutrino flavor (default is nutau): 12 for nue 14 for numu 16 for nutau'
                       )

    # Energy arguments
    parser.add_argument('-e', 
                        dest='energy',
                        default='',
                        help='Energy in GeV if numerical value greater than 0 passed.\n\
                              Sprectral index if numerical value less than or equal to 0 passed.\n\
                              Else path to CDF to sample from.'
                       )
    parser.add_argument('--e_min',
                        type=float,
                        default=1e6,
                        help='Minimum energy to sample from if doing powerlaw sampling'
                       )
    parser.add_argument('--e_max',
                        type=float,
                        default=1e9,
                        help='Minimum energy to sample from if doing powerlaw sampling'
                       )
    
    # Angular arguments
    parser.add_argument('-t', 
                        dest='theta',
                        default='',
                        help='nadir angle in degrees if numerical value(0 is through the core).\n\
                              "Range" if you want to sample from a range of thetas. Use --th_min and --th_max to specify lims.'
                       )
    parser.add_argument('--th_max',
                        dest='th_max',
                        type=float, 
                        default=90,
                        help='If doing a theta range, maximum theta value. Default 90, i.e. skimming'
                       )
    parser.add_argument('--th_min', 
                        type=float,
                        default=0,
                        help='If doing a theta range, minimum theta value. Default 0, i.e. through the core'
                       )
    
    # Saving arguments
    parser.add_argument('--savedir', 
                        dest='savedir', 
                        type=str, 
                        default='', 
                        help="If saving output, provide a path here, if not specified, output will be printed"
                       )
    parser.add_argument('--prefix',
                        dest='prefix', 
                        default='',
                        type=str,
                        help='(Optional) argument to give specify name of outdir'
                       )

    # cross section args
    parser.add_argument('--xs',
                        dest='xs_model',
                        type=str,
                        default='dipole',
                        help="Enter 'CSMS' if you would like to run the simulation with a pQCD xs model"
                       )

    # Body specification args
    parser.add_argument('--body',
                        dest='body',
                        default='earth',
                        help="body to use"
                       )
    parser.add_argument('--water',
                        dest='water',
                        type=float, default=0,
                        help="If you you would like to add a water layer to the Earth model, enter it here in km."
                       )
    parser.add_argument('--radius',
                        dest='radius',
                        type=float,
                        default=0,
                        help="Raise this flag if you want to turn off tau losses. In this case, taus will decay at rest."
                       )
    parser.add_argument('--depth',
                        dest='depth',
                        type=float,
                        default=0,
                        help="Depth of the detector in km."
                       )
    
    # Options
    parser.add_argument('--no_losses',
                        dest='no_losses', 
                        default=True,
                        action='store_false',
                        help="Raise this flag if you want to turn off tau losses. In this case, taus will decay at rest."
                       )
    parser.add_argument('--no_secondaries',
                        default=False,
                        action='store_true',
                        help="Raise this flag to turn off secondaries"
                       )
    parser.add_argument('-d', 
                        dest='debug', 
                        default=False, 
                        action='store_true', 
                        help='Do you want to print out debug statments? If so, raise this flag'
                       ) 

    args = parser.parse_args()
    return args

def run_MC(eini, thetas, body, xs, tracks, TR_specs):
    r'''
    Main simulation code. Propagates a flux of neutrinos and returns or
    saves the outgoing particles


    Returns:
        np.recarray of the outgoing leptons
    '''

    propagated_stack     = []
    inds_left            = list(range(TR_specs['nevents']))
    output               = []
    counter              = 0
    iter_energies        = list(eini)
    iter_positions       = list(np.zeros(TR_specs['nevents']))
    iter_particleID      = np.ones(TR_specs['nevents'], dtype=int)*TR_specs['flavor']
    iter_ChargedPosition = list(np.zeros(TR_specs['nevents']))
    iter_nCC             = list(np.zeros(TR_specs['nevents']))
    iter_nNC             = list(np.zeros(TR_specs['nevents']))
    rand                 = TR_specs['rand']

    # Run the algorithm
    # All neutrinos are propagated until either exiting or undergoing a CC interaction.
    # All CC interactions are handled together, and then the next iteration occurs
    # This repeats until all leptons have reached the total distance
    t0 = time.time()

    while inds_left:
        counter += 1
        cc_stack = []

        for j in range(len(inds_left) - 1, -1, -1):
            i = inds_left[j] #Unique event index

            particle = Particle(iter_particleID[i], iter_energies[i], 
                                thetas[i], iter_positions[i], i, rand.randint(low=1e9),
                                iter_ChargedPosition[i], xs, not TR_specs['no_secondaries'])
            my_track = tracks[thetas[i]]
            out = Propagate(particle, my_track, body)

            iter_nCC[i]+=out.nCC
            iter_nNC[i]+=out.nNC

            if (out.survived==False):
                #these muons were absorbed. we record them in the output with outgoing energy 0
                output.append((eini[ind], 0., thetas[ind], inds_left[j], iter_nCC[ind], iter_nNC[ind], out.ID))
                iter_positions[int(out.index)] = float(out.position)
                del inds_left[j]
                del out
            elif (out.isCC):
                current_distance=my_track.x_to_d(out.position)*body.radius
                current_x = out.position
                total_distance=my_track.x_to_d(1.)*body.radius
                #current_density=body.get_density(my_track.x_to_r(out.position))
                current_density=body.get_average_density(my_track.x_to_r(out.position))
                cc_stack.append((float(out.energy), current_x, float(current_distance), int(out.index),
                 int(out.ID), 0, float(total_distance), float(current_density), int(out.ID)))
                del out
            else:
                ind = int(out.index)
                if ind != i:
                    message += "Index mismatch: {} {}".format(ind, i)
                    raise RuntimeError('Index mismatch -- particles are getting jumbled somewhere (thats a bad thing)')
                output.append((eini[ind], float(out.energy), thetas[ind], inds_left[j], iter_nCC[ind], iter_nNC[ind], out.ID))
                iter_positions[int(out.index)] = float(out.position)
                if not TR_specs['no_secondaries']:
                    basket = out.basket
                    for sec in basket:
                        sec_particle = Particle(sec['ID'], sec['energy'], thetas[ind], sec['position'], inds_left[j], rand.randint(low=1e9),
                                                            0.0, xs=xs,secondaries=False)
                        sec_out      = Propagate(sec_particle, my_track, body)
                        if(sec_out.isCC):
                            output.append((sec_out.energy, 0.0, thetas[ind], inds_left[j], sec_out.nCC, sec_out.nNC, -sec_out.ID))
                        else:
                            output.append((sec_out.initial_energy, sec_out.energy, thetas[ind], 
                                               inds_left[j], sec_out.nCC, sec_out.nNC, -sec_out.ID))
                        del sec_particle
                    del out.basket
                    del basket
                del inds_left[j]
                del out           

        if (len(cc_stack) > 0):
            #if debug:
            #    message += "{} events passed to MMC in loop iteration {}\n".format(len(cc_stack), counter)
            EventCollection = DoAllCCThings(cc_stack, xs, not TR_specs['no_losses'])
            for event in EventCollection:
                iter_positions[int(event[3])] = float(event[1])
                iter_energies[int(event[3])] = float(event[0])
                iter_particleID[int(event[3])] = int(event[4])
                iter_ChargedPosition[int(event[3])] = float(event[5])
                del event
    #if debug:
    #    print("Simulating {} events at {} degrees took {} seconds.".format(nevents, theta, time.time() - t0))

    output = np.array(output, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('CDF_index', float), ('nCC', int), ('nNC', int), ('PDG_Encoding', int)])
    output['Theta'] *= 180. / np.pi #Give theta in degrees to user

    return output

if __name__ == "__main__": # pragma: no cover

    args = initialize_parser()
    TR_specs = {}
    TR_specs['seed']           = args.seed
    TR_specs['nevents']        = args.nevents
    TR_specs['flavor']         = args.flavor
    TR_specs['energy']         = args.energy
    TR_specs['e_min']          = args.e_min
    TR_specs['e_max']          = args.e_max
    TR_specs['theta']          = args.theta
    TR_specs['th_min']         = args.th_min
    TR_specs['th_max']         = args.th_max
    TR_specs['base_savedir']   = args.savedir
    TR_specs['water']          = args.water
    TR_specs['xs_model']       = args.xs_model
    TR_specs['no_losses']      = args.no_losses
    TR_specs['body']           = args.body
    TR_specs['radius']         = args.radius
    TR_specs['depth']          = args.depth
    TR_specs['no_secondaries'] = args.no_secondaries
    TR_specs['prefix']         = args.prefix
    TR_specs['debug']          = ''
    
    if TR_specs['nevents']<=0:
        raise ValueError("We need to simulate at least one event, c'mon y'all")

    # Set up outdir and make seed if not provided
    if TR_specs['base_savedir']:
        if not os.path.isdir(TR_specs['base_savedir']):
            raise ValueError('Savedir %s does not exist' % TR_specs['base_savedir'])
        from taurunner.modules import setup_outdir
        TR_specs = setup_outdir(TR_specs)
    
    # Set up a random state
    rand = np.random.RandomState(TR_specs['seed'])
    TR_specs['rand'] = rand

    # Make an array of injected energies
    from taurunner.modules import make_initial_e
    eini = make_initial_e(TR_specs, rand=rand)

    # Make an array of injected incident angles
    from taurunner.modules import make_initial_thetas
    thetas = make_initial_thetas(TR_specs, rand=rand)

    # Make body
    # TODO Make this not suck. i.e. make construct body more comprehensive
    from taurunner.modules import construct_body
    body = construct_body(TR_specs)

    # Premake all necessary tracks in case of redundancies
    # TODO Make it so that you can pass radial tracks too
    tracks  = {theta:Chord(theta=theta, depth=TR_specs['depth']/body.radius) for theta in set(thetas)}
    # TODO Make cross section obect
    xs = CrossSections(TR_specs['xs_model'])

    result = run_MC(eini, thetas, body, xs, tracks, TR_specs)

    if TR_specs['base_savedir']:
        base_fname = '%s/%s/%s' % (TR_specs['base_savedir'], TR_specs['prefix'], TR_specs['prefix'])
        np.save(base_fname+'.npy', result)
        with open(base_fname + '.json', 'w') as f:
            TR_specs.pop('rand')
            json.dump(TR_specs[mm], f)
        if TR_specs['debug']:
            print(message)
    else:
        if TR_specs['debug']:
            print(message)
        try:
            from tabulate import tabulate
            headers = list(result.dtype.names)
            out_table = tabulate(result, headers, tablefmt="fancy_grid")
            print(out_table)
        except ImportError:
            print("Outgoing Particles: ")
            print(result)
