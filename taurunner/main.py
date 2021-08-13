#!/usr/bin/env python

import os, sys, json
import proposal as pp

from taurunner.utils import units, sample_powerlaw, is_floatable, make_propagator
import taurunner.track as track
from taurunner.body import *
from taurunner.cross_sections import CrossSections
from taurunner.Casino import *
from taurunner.particle import Particle

def initialize_parser(): # pragma: no cover
    import argparse
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
    parser.add_argument('--track',
                        dest='track',
                        type=str,
                        default='Chord',
                        help='Track type to use. curently only radial and chord trajectories supported.'
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
    parser.add_argument('--save', 
                        dest='save', 
                        type=str, 
                        default='', 
                        help="If saving output, provide a path here, if not specified, output will be printed"
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
                        default=False,
                        action='store_true',
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

def run_MC(eini: np.ndarray, 
           thetas: np.ndarray,
           body: Body, 
           xs: CrossSections, 
           tracks: dict, 
           TR_specs: dict, 
           propagator: pp.Propagator,
           rand: np.random.RandomState = None
          ) -> np.ndarray:
    r'''
    Main simulation code. Propagates an ensemble of initial states and returns the output
    Params
    ______
    eini       : array containing the initial energies of particles to simulate
    thetas     : array containing the incoming angles of particles to simulate
    body       : taurunner Body object in which to propagate the particles
    xs         : taurunner CrossSections object for interactions
    tracks     : dictionary whose keys are angles and who values are taurunner Track objects
    TR_specs   : dictionary specifiying additional TR options
    propagator : PROPOSAL propagator object for charged lepton propagation
    Returns
    _______
    output : Array containing the output information of the MC. This includes 
             initial and final energies, incident incoming angles, number of CC and NC interactions,
             and particle type (PDG convention)
             
    '''
    output      = []
    energies    = list(eini)
    particleIDs = np.ones(TR_specs['nevents'], dtype=int)*TR_specs['flavor']

    if rand is None:
        rand = np.random.RandomState()
    secondary_basket = []
    idxx             = []

    # Run the algorithm
    # All neutrinos are propagated until exiting as tau neutrino or taus.
    # If secondaries are on, then each event has a corresponding secondaries basket
    # which are propagated all at once in the end.
    for i in range(TR_specs['nevents']):
        particle = Particle(particleIDs[i], energies[i], thetas[i], 0.0, rand, xs,
                            propagator, not TR_specs['no_secondaries'], TR_specs['no_losses'])
        
        my_track = tracks[thetas[i]]
        out = Propagate(particle, my_track, body)
    
        if (out.survived==False):
            #these muons/electrons were absorbed. we record them in the output with outgoing energy 0
            output.append((energies[i], 0., thetas[i], out.nCC, out.nNC, out.ID))
        else:
            output.append((energies[i], float(out.energy), thetas[i], out.nCC, out.nNC, out.ID))
        if not TR_specs['no_secondaries']:
            secondary_basket.append(np.asarray(out.basket))
            idxx = np.hstack([idxx, [i for _ in out.basket]])
            del out.basket
        del out     
        del particle
    idxx = idxx.astype(np.int32)
    if not TR_specs['no_secondaries']:    
        #make muon propagator
        secondary_basket = np.concatenate(secondary_basket)
        #ids = np.unique([s['ID'] for s in secondary_basket])
        sec_prop = {ID:make_propagator(ID, body, TR_specs['xs_model']) for ID in [-12, -14]}
        for sec, i in zip(secondary_basket, idxx):
            sec_particle = Particle(sec['ID'], sec['energy'], thetas[i], sec['position'], rand,
                                    xs=xs, proposal_propagator=sec_prop[sec['ID']], secondaries=False, no_losses=False)
            sec_out      = Propagate(sec_particle, my_track, body)
            if(not sec_out.survived):
                output.append((sec_out.energy, 0.0, thetas[i], sec_out.nCC, sec_out.nNC, sec_out.ID))
            else:
                output.append((sec_out.initial_energy, sec_out.energy, thetas[i], 
    				                 sec_out.nCC, sec_out.nNC, sec_out.ID))
            del sec_particle
            del sec_out     
        
    output = np.array(output, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('nCC', int), ('nNC', int), ('PDG_Encoding', int)])
    output['Theta'] = np.degrees(output['Theta'])
                       
    return output

if __name__ == "__main__": # pragma: no cover

    args = initialize_parser()
    TR_specs = {}
    TR_specs['seed']           = args.seed
    TR_specs['nevents']        = args.nevents
    TR_specs['flavor']         = args.flavor
    TR_specs['track']          = args.track
    TR_specs['energy']         = args.energy
    TR_specs['e_min']          = args.e_min
    TR_specs['e_max']          = args.e_max
    TR_specs['theta']          = args.theta
    TR_specs['th_min']         = args.th_min
    TR_specs['th_max']         = args.th_max
    TR_specs['save']           = args.save
    TR_specs['water']          = args.water
    TR_specs['xs_model']       = args.xs_model
    TR_specs['no_losses']      = args.no_losses
    TR_specs['body']           = args.body
    TR_specs['radius']         = args.radius
    TR_specs['depth']          = args.depth
    TR_specs['no_secondaries'] = args.no_secondaries
    TR_specs['debug']          = ''

    if TR_specs['nevents']<=0:
        raise ValueError("We need to simulate at least one event, c'mon y'all")

    # Set up outdir and make seed if not provided
    if TR_specs['save']:
        savedir = '/'.join(TR_specs['save'].split('/')[:-1])
        if not os.path.isdir(savedir):
            raise ValueError('Savedir %s does not exist' % TR_specs['save'])
    
    # Set up a random state
    rand = np.random.RandomState(TR_specs['seed'])

    # Make body
    # TODO Make this not suck. i.e. make construct body more comprehensive
    from taurunner.utils import construct_body
    body = construct_body(TR_specs)

    # Make an array of injected energies
    from taurunner.utils.make_initial_e import make_initial_e
    eini = make_initial_e(TR_specs, rand=rand)

    # Make an array of injected incident angles
    from taurunner.utils.make_initial_thetas import make_initial_thetas
    thetas = make_initial_thetas(TR_specs, rand=rand)

    from taurunner.utils.make_tracks import make_tracks
    tracks = make_tracks(TR_specs, thetas, body.radius)

    xs = CrossSections(TR_specs['xs_model'])

    prop = make_propagator(TR_specs['flavor'], body, TR_specs['xs_model'])

    result = run_MC(eini, thetas, body, xs, tracks, TR_specs, prop, rand=rand)

    if TR_specs['save']:
        if '.npy' not in TR_specs['save']:
            TR_specs['save'] += '.npy'
        np.save(TR_specs['save'], result)
        with open(TR_specs['save'].replace('npy', 'json'), 'w') as f:
            json.dump(TR_specs, f)
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
