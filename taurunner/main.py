#!/usr/bin/env python

import numpy as np
import os
import json

from typing import Union

from Casino import run_MC, CrossSections, XSModel

def string_or_float(v: str) -> Union[str, float]:
    """
    Return float version of arg if floatable, else return arg
    """
    try:
        return float(v)
    except ValueError:
        return v

def posint(v: str) -> int:
    """
    return integer version of str assuring it is positive
    """
    i = int(v)
    if i<=0:
        raise ValueError("We need to simulate at least one event, c'mon y'all")
    return i

def initialize_parser(): # pragma: no cover
    r'''
    Helper function to parse command-line arguments
    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        dest='seed',
        type=int,
        default=None,
        help='just an integer seed to help with output file names'
    )
    parser.add_argument(
        '-n',
        dest='nevents',
        type=posint,
        required=True,
        help='how many events do you want?'
    )
    parser.add_argument(
        '--flavor',
        dest='flavor',
        type=int, default=16,
        help='neutrino flavor (default is nutau): 12 for nue 14 for numu 16 for nutau'
    )
    parser.add_argument(
        '--track',
        dest='track',
        type=str,
        default='chord',
        help='Track type to use. curently only radial and chord trajectories supported.'
    )
    # Energy arguments
    parser.add_argument(
        '-e',
        dest='energy',
        type=string_or_float,
        required=True,
        help='Energy in GeV if numerical value greater than 0 passed.\n\
             Sprectral index if numerical value less than or equal to 0 passed.\n\
             Else path to CDF to sample from.'
    )
    parser.add_argument(
        '--e_min',
        type=float,
        default=1e15,
        help='Minimum energy to sample from if doing powerlaw sampling (eV)'
    )
    parser.add_argument(
        '--e_max',
        type=float,
        default=1e18,
        help='Maximum energy to sample from if doing powerlaw sampling (eV)'
    )
    
    # Angular arguments
    parser.add_argument(
        '-t',
        dest='theta',
        type=string_or_float,
        required=True,
        help='nadir angle in degrees if numerical value(0 is through the core).\n\
            "range" if you want to sample from a range of thetas. Use --th_min and --th_max to specify lims.'
    )
    parser.add_argument(
        '--th_max',
        dest='th_max',
        type=float,
        default=90,
        help='If doing a theta range, maximum theta value. Default 90, i.e. skimming'
    )
    parser.add_argument(
        '--th_min',
        type=float,
        default=0,
        help='If doing a theta range, minimum theta value. Default 0, i.e. through the core'
    )

    # Saving arguments
    parser.add_argument(
        '--save',
        dest='save',
        type=str,
        default='',
        help="If saving output, provide a path here, if not specified, output will be printed"
    )

    # cross section args
    parser.add_argument(
        '--xs',
        dest='xs_model',
        type=str,
        default='dipole',
        help="Enter 'CSMS' if you would like to run the simulation with a pQCD xs model"
    )

    # Body specification args
    parser.add_argument(
        '--body',
        dest='body',
        default='earth',
        help="body to use"
    )
    parser.add_argument(
        '--water',
        dest='water',
        type=float, default=0,
        help="If you you would like to add a water layer to the Earth model, enter it here in km."
    )
    parser.add_argument(
        '--radius',
        dest='radius',
        type=float,
        default=0,
        help="Radius of some custom body"
    )
    parser.add_argument(
        '--depth',
        dest='depth',
        type=float,
        default=0,
        help="Depth of the detector in km."
    )
    # Options
    parser.add_argument(
        '--no_losses',
        dest='no_losses',
        default=False,
        action='store_true',
        help="Raise this flag if you want to turn off tau losses. In this case, taus will decay at rest."
    )
    parser.add_argument(
        '--no_secondaries',
        default=False,
        action='store_true',
        help="Raise this flag to turn off secondaries"
    )
    parser.add_argument(
        '--e_cut',
        dest='e_cut',
        default=0.0,
        type=float,
        help='Energy at which to stop the propagation.'
    )

    args = parser.parse_args()
    return args

def main() -> None:
    args = initialize_parser()

    np.random.seed(args.seed)

    # Construct the body
    from taurunner.utils import construct_body, make_initial_e, make_initial_thetas
    body = construct_body(args.body, args.water, radius=args.radius)

    eini = make_initial_e(
        args.nevents,
        args.energy,
        e_min=args.e_min,
        e_max=args.e_max,
    )

    if(args.theta=='range'):
        theta = (args.th_min, args.th_max)
    else:
        theta = args.theta
    thetas = make_initial_thetas(args.nevents, theta, track_type=args.track)

    # Set up the cross sections
    xs = CrossSections(getattr(XSModel, args.xs_model.upper()))

    condition = None
    if args.e_cut:
        condition = lambda particle: (particle.energy <= args.e_cut*units.GeV and abs(particle.ID) in [12,14,16])

    result = run_MC(
        eini,
        thetas,
        body,
        xs,
        no_secondaries=args.no_secondaries,
        no_losses=args.no_losses,
        flavor=args.flavor,
        condition=condition,
        depth=args.depth,
        track_type=args.track
    )

    if args.save:
        savedir = '/'.join(args.save.split('/')[:-1])
        if not os.path.isdir(savedir):
            raise ValueError('Savedir %s does not exist' % args.save)
        if '.npy' not in args.save:
            args.save += '.npy'
        np.save(args.save, result)
        with open(args.save.replace('npy', 'json'), 'w') as f:
            json.dump(vars(args), f)
    else:
        from tabulate import tabulate
        headers = list(result.dtype.names)
        out_table = tabulate(result, headers, tablefmt="fancy_grid")
        print(out_table)

if __name__ == "__main__": # pragma: no cover
    main()
