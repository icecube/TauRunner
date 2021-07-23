#!/usr/bin/env python

import os, sys
import json
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import argparse

from taurunner.modules import units, cleanup_outdir, sample_powerlaw, is_floatable
#from taurunner.modules import units, make_outdir, todaystr, cleanup_outdir
from taurunner.track import Chord
from taurunner.body import *
from taurunner.cross_sections import CrossSections
from taurunner.Casino import *
import proposal as pp

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
                        default=False,
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

def run_MC(eini, thetas, body, xs, tracks, TR_specs, propagator):
    r'''
    Main simulation code. Propagates a flux of neutrinos and returns or
    saves the outgoing particles


    Returns:
        np.recarray of the outgoing leptons
    '''
    output                 = []
    energies               = list(eini)
    particleIDs            = np.ones(TR_specs['nevents'], dtype=int)*TR_specs['flavor']
    rand                   = TR_specs['rand']
    proposal_lep           = pp.particle.DynamicData(pp.particle.TauMinusDef().particle_type)
    proposal_lep.position  = pp.Vector3D(0, 0, 0)
    proposal_lep.direction = pp.Vector3D(0, 0, 1)

    # Run the algorithm
    # All neutrinos are propagated until exiting as tau neutrino or taus.
    # If secondaries are on, then each event has a corresponding secondaries basket
    # which are propagated all at once in the end.
    t0 = time.time()
    for i in range(TR_specs['nevents']):
        particle = Particle(particleIDs[i], energies[i], thetas[i], 0.0, rand.randint(low=1e9),
    			    xs, propagator, proposal_lep, not TR_specs['no_secondaries'], TR_specs['no_losses'])
        
        my_track = tracks[thetas[i]]
        out = Propagate(particle, my_track, body)
    
        if (out.survived==False):
            #these muons were absorbed. we record them in the output with outgoing energy 0
            output.append((energies[i], 0., thetas[i], out.nCC, out.nNC, out.ID))
            del out
        else:
            output.append((energies[i], float(out.energy), thetas[i], out.nCC, out.nNC, out.ID))
        if not TR_specs['no_secondaries']:
            basket = out.basket
            for sec in basket:
                sec_particle = Particle(sec['ID'], sec['energy'], thetas[i], sec['position'], rand.randint(low=1e9),
                                        xs=xs, proposal_propagator=None, proposal_lep=None, secondaries=False, no_losses=True)
                sec_out      = Propagate(sec_particle, my_track, body)
                if(not sec_out.survived):
                    output.append((sec_out.energy, 0.0, thetas[i], sec_out.nCC, sec_out.nNC, -sec_out.ID))
                else:
                    output.append((sec_out.initial_energy, sec_out.energy, thetas[i], 
    				                 sec_out.nCC, sec_out.nNC, -sec_out.ID))
                del sec_particle
            del out.basket
            del basket
        del out           
    output = np.array(output, dtype = [('Eini', float), ('Eout',float), ('Theta', float), ('nCC', int), ('nNC', int), ('PDG_Encoding', int)])
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
    rand = np.random.RandomState(TR_specs['seed'])
    TR_specs['rand'] = rand
    if TR_specs['nevents']<=0:
        raise ValueError("We need to simulate at least one event, c'mon y'all")

    # Set up outdir and make seed if not provided
    if TR_specs['base_savedir']:
        if not os.path.isdir(TR_specs['base_savedir']):
            raise ValueError('Savedir %s does not exist' % TR_specs['base_savedir'])
        from taurunner.modules import setup_outdir
        TR_specs = setup_outdir(TR_specs)

    try:        
        # Make injected energies
        if is_floatable(TR_specs['energy']):
            e = float(TR_specs['energy'])
            if e<=0:
                eini = sample_powerlaw(rand, TR_specs['e_min'], TR_specs['e_max'], e + 1, size=TR_specs['nevents'])*units.GeV
            else:
                eini = np.full(TR_specs['nevents'], e)*units.GeV
        else:
            if not os.path.isfile(TR_specs['energy']):
                raise RuntimeError("GZK CDF Spline file does not exist")
            # sample initial energies and incoming angles from GZK parameterization
            cdf_indices = rand.uniform(low=0., high=1.,size=TR_specs['nevents'])
            cdf         = np.load(TR_specs['energy'], allow_pickle=True).item()
            eini        = cdf(cdf_indices)*units.GeV
            if TR_specs['debug']:
                message+="Sampled {} events from the GZK flux\n".format(nevents)

        # Make injected thetas
        if is_floatable(TR_specs['theta']):
            t = float(TR_specs['theta'])
            if t<0 or t>90:
                raise ValueError('Angles must be between 0 and 90')
            thetas = np.radians(np.ones(TR_specs['nevents'])*t)
        else:
            if TR_specs['theta']=='range':
                if TR_specs['th_min']<0 or TR_specs['th_max']>90:
                    raise ValueError('Angles must be between 0 and 90')
                costh_min = np.cos(np.radians(TR_specs['th_max']))
                costh_max = np.cos(np.radians(TR_specs['th_min']))
                thetas = np.arccos(rand.uniform(costh_min, costh_max, TR_specs['nevents']))
            else:
                raise ValueError('theta sampling %s not suppoorted' % TR_specs['theta'])

        # Make body
        # TODO Make this not suck. i.e. make construct body more comprehensive
        from taurunner.modules import construct_body
        body = construct_body(TR_specs)

        # Premake all necessary tracks in case of redundancies
        # Make it so that you can pass radial tracks too
        tracks  = {theta:Chord(theta=theta, depth=TR_specs['depth']/body.radius) for theta in set(thetas)}
        rr = np.linspace(0, 1, 50)
        #print(np.asarray([body.get_average_density(r) for r in rr])/units.gr*units.cm**3)
        # Make cross section obect
        xs = CrossSections(TR_specs['xs_model'])

        #Make proposal object
        #define geometry
        sec_def = pp.SectorDefinition()
        sec_def.medium = pp.medium.Ice(1.0)
        sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), 1e20, 0)
        sec_def.particle_location = pp.ParticleLocation.inside_detector        
        sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere
        sec_def.crosssection_defs.brems_def.lpm_effect = True
        sec_def.crosssection_defs.epair_def.lpm_effect = True
        
        sec_def.cut_settings.ecut = 500
        sec_def.cut_settings.vcut = 0.1
        
        sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.BlockDurandHa
        
        interpolation_def = pp.InterpolationDef()
        interpolation_def.path_to_tables = "/home/isafa/.local/share/PROPOSAL/tables"
        interpolation_def.path_to_tables_readonly = "/home/isafa/.local/share/PROPOSAL/tables"
        interpolation_def.nodes_cross_section = 200
        
        #define propagator
        prop = pp.Propagator(
                particle_def=pp.particle.TauMinusDef(),
                sector_defs=[sec_def],
                detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
                interpolation_def=interpolation_def
        )

        result = run_MC(eini, thetas, body, xs, tracks, TR_specs, prop)

        if TR_specs['base_savedir']:
            base_fname = '%s/%s/%s' % (TR_specs['base_savedir'], TR_specs['prefix'], TR_specs['prefix'])
            np.save(base_fname+'.npy', result)
            with open(base_fname + '.json', 'w') as f:
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
    
    except Exception as e:
        #if TR_specs['base_savedir']:
        #    os.rmdir('%s/%s' % (TR_specs['base_savedir'], TR_specs['prefix']))
        raise e
