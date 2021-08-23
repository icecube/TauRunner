import os
import numpy as np
import proposal as pp
from scipy.optimize import ridder
from importlib.resources import path
from scipy.integrate import quad
from taurunner.utils import units

ID_2_name = {16: 'TauMinusDef', -16: 'TauPlusDef',
             14: 'MuMinusDef',  -14: 'MuPlusDef',
             12: 'EMinusDef',   -12: 'EPlusDef',}

def segment_body(body, granularity=0.5):
    descs = []
    for xi, xf in zip(body.layer_boundaries[:-1], body.layer_boundaries[1:]):
        gran = np.abs(body.get_density(xf, right=True) - body.get_density(xi, right=False)) / body.get_density(xi, right=False)
        if gran<granularity:  # the percent difference within a layer is small
            descs.append((xi, xf, body.get_average_density(0.5*(xi+xf))))
        else:
            end   = 0.0
            start = xi
            while end<xf:
                s_density = body.get_density(start, right=False)
                func      = lambda x: np.abs(body.get_density(x, right=True)-s_density)/s_density-granularity
                end       = ridder(func, xi, xf)
                if end>xf: # you've left the layer so only go to end of layer
                    end = xf
                wrap = lambda x: body.get_density(x, right=True)
                I = quad(wrap, start, end, full_output=1)
                avg_density  = I[0]/(end-start)
                descs.append((start, end, avg_density))
                start = end
    return descs

def make_propagator(ID, body, xs_model='dipole', granularity=0.5):
    
    if(ID in [12, -12]):
        return None
    
    particle_def              = getattr(pp.particle, ID_2_name[ID])()

    #define how many layers of constant density we need for the tau
    descs = segment_body(body, granularity)
    #make the sectors
    sec_defs = [make_sector(d/units.gr*units.cm**3, s*body.radius/units.meter, e*body.radius/units.meter, xs_model) for s, e, d in descs]

    with path('taurunner.resources.proposal_tables', 'tables.txt') as p:
        tables_path = str(p).split('tables.txt')[0]
    
    #define interpolator
    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = tables_path
    interpolation_def.path_to_tables_readonly = tables_path
    interpolation_def.nodes_cross_section = 199

    #define propagator -- takes a particle definition - sector - detector - interpolator
    prop = pp.Propagator(particle_def=particle_def,
                         sector_defs=sec_defs,
                         detector=pp.geometry.Sphere(pp.Vector3D(), body.radius/units.cm, 0),
                         interpolation_def=interpolation_def)

    return prop


def make_sector(density, start, end, xs_model):
    sec_def = pp.SectorDefinition()
    sec_def.medium = pp.medium.Ice(density)
    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), end, start)
    sec_def.particle_location = pp.ParticleLocation.inside_detector        
    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere
    sec_def.crosssection_defs.brems_def.lpm_effect = True
    sec_def.crosssection_defs.epair_def.lpm_effect = True
    
    sec_def.cut_settings.ecut = 1e4*1e3
    sec_def.cut_settings.vcut = 1e-3
    sec_def.do_continuous_randomization = False
    if(xs_model=='dipole'):
        sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.BlockDurandHa
    else:
        sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.AbramowiczLevinLevyMaor97

    return sec_def