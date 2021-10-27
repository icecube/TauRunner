#!/usr/bin/env python

import os, sys, json
import proposal as pp

from taurunner.utils import units, make_propagator
from taurunner import track
from taurunner.body import *
from taurunner.cross_sections import CrossSections
from taurunner.Casino import *
from taurunner.particle import Particle
from taurunner.utils.make_track import make_track

def initialize_parser(): # pragma: no cover
    r'''
    Helper function to parse command-line arguments
    '''
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
    parser.add_argument('--flavor', 
                        dest='flavor', 
                        type=int, default=16,
                        help='neutrino flavor (default is nutau): 12 for nue 14 for numu 16 for nutau'
                       )
    parser.add_argument('--track',
                        dest='track',
                        type=str,
                        default='chord',
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
                        default=1e15,
                        help='Minimum energy to sample from if doing powerlaw sampling (eV)'
                       )
    parser.add_argument('--e_max',
                        type=float,
                        default=1e18,
                        help='Maximum energy to sample from if doing powerlaw sampling (eV)'
                       )
    
    # Angular arguments
    parser.add_argument('-t', 
                        dest='theta',
                        default='',
                        help='nadir angle in degrees if numerical value(0 is through the core).\n\
                              "range" if you want to sample from a range of thetas. Use --th_min and --th_max to specify lims.'
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
    parser.add_argument('--e_cut',
                        dest='e_cut',
                        default=0.0,
                        type=float,
                        help='Energy at which to stop the propagation.'
                       )

    args = parser.parse_args()
    return args

def run_MC(eini: np.ndarray, 
           thetas: np.ndarray,
           body: Body, 
           xs: CrossSections, 
           propagator: pp.Propagator,
           rand: np.random.RandomState = None,
           no_secondaries: bool = False,
           flavor: int = 16,
           no_losses: bool = False,
           condition=None,
           depth=0.0,
           track_type='chord'
          ) -> np.ndarray:
    r'''
    Main simulation code. Propagates an ensemble of initial states and returns the output
    Params
    ______
    eini       : array containing the initial energies of particles to simulate
    thetas     : array containing the incoming angles of particles to simulate
    body       : taurunner Body object in which to propagate the particles
    xs         : taurunner CrossSections object for interactions
    propagator : PROPOSAL propagator object for charged lepton propagation
    Returns
    _______
    output : Array containing the output information of the MC. This includes 
             initial and final energies, incident incoming angles, number of CC and NC interactions,
             and particle type (PDG convention)
             
    '''
    nevents     = len(eini)
    output      = []
    energies    = list(eini)
    particleIDs = np.ones(nevents, dtype=int)*flavor
    xs_model    = xs.model
    prev_th     = thetas[0]
    prev_track  = make_track(prev_th)
    if rand is None:
        rand = np.random.RandomState()
    secondary_basket = []
    idxx             = []    
    my_track  = None
    prv_theta = np.nan
    # Run the algorithm

    # All neutrinos are propagated until exiting as tau neutrino or taus.
    # If secondaries are on, then each event has a corresponding secondaries basket
    # which are propagated all at once in the end.
    for i in range(nevents):
        cur_theta = thetas[i]
        cur_e     = energies[i]

        if (cur_theta!=prv_theta and track_type=='chord') or my_track is None: # We need to make a new track
            my_track = getattr(track, track_type)(theta=cur_theta, depth=depth)
        particle = Particle(particleIDs[i], 
                            cur_e, 
                            0.0 ,
                            rand, 
                            xs,
                            propagator, 
                            not no_secondaries, 
                            no_losses
        )
        
        out      = Propagate(particle, my_track, body, condition=condition)
    
        if (out.survived==False):
            #this muon/electron was absorbed. we record it in the output with outgoing energy 0
            output.append((cur_e, 0., cur_theta, out.nCC, out.nNC, out.ID, i, out.position))
        else:
            #this particle escaped
            output.append((cur_e, float(out.energy), cur_theta, out.nCC, out.nNC, out.ID, i, out.position))
        if not no_secondaries:
            #store secondaries to propagate later
            secondary_basket.append(np.asarray(out.basket))
            #keep track of parent
            idxx = np.hstack([idxx, [i for _ in out.basket]])
        prv_theta = cur_theta
        del out
        del particle
    idxx = np.asarray(idxx).astype(np.int32)
    if not no_secondaries:    
        #make muon propagator
        secondary_basket = np.concatenate(secondary_basket)
        sec_prop         = {ID:make_propagator(ID, body, xs_model) for ID in [-12, -14]}
        for sec, i in zip(secondary_basket, idxx):
            cur_theta = thetas[i]
            if cur_theta!=prv_theta and track_type=='chord': # We need to make a new track
                my_track = getattr(track, 'chord')(theta=cur_theta, depth=depth)
            sec_particle = Particle(
                                    sec['ID'], 
                                    sec['energy'],
                                    sec['position'],
                                    rand,
                                    xs=xs, 
                                    proposal_propagator=sec_prop[sec['ID']], 
                                    secondaries=False, 
                                    no_losses=False
                                   )
            sec_out = Propagate(sec_particle, my_track, body, condition=condition)
            if(not sec_out.survived):
                output.append((sec_out.initial_energy, 0.0, cur_theta, sec_out.nCC, sec_out.nNC, sec_out.ID, i, sec_out.position))
            else:
                output.append((sec_out.initial_energy, sec_out.energy, cur_theta, 
    	    		                 sec_out.nCC, sec_out.nNC, sec_out.ID, i, sec_out.position))
            prv_theta = cur_theta
            del sec_particle
            del sec_out
    output = np.array(output, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('nCC', int), ('nNC', int), ('PDG_Encoding', int), ('event_ID', int), ('final_position', float)])
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
    eini = make_initial_e(TR_specs['nevents'], TR_specs['energy'],
                          e_min=TR_specs['e_min'], e_max=TR_specs['e_max'], rand=rand)

    # Make an array of injected incident angles
    from taurunner.utils.make_initial_thetas import make_initial_thetas
    if(TR_specs['theta']=='range'):
        theta = (TR_specs['th_min'], TR_specs['th_max'])
    else:
        theta = TR_specs['theta']
    thetas = make_initial_thetas(TR_specs['nevents'], theta, rand=rand, track_type=TR_specs['track'])
    sorter = np.argsort(thetas)
    thetas = thetas[sorter]
    
    xs = CrossSections(TR_specs['xs_model'])

    prop = make_propagator(TR_specs['flavor'], body, TR_specs['xs_model'])
    if args.e_cut:
        condition = lambda particle: (particle.energy <= args.e_cut*units.GeV and abs(particle.ID) in [12,14,16])
    else:
        condition = None

    result = run_MC(eini, thetas, body, xs, prop, rand=rand,
                    no_secondaries=TR_specs['no_secondaries'], no_losses=TR_specs['no_losses'],
                    flavor=TR_specs['flavor'], condition=condition, depth=TR_specs['depth'], track_type=TR_specs['track'])

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
