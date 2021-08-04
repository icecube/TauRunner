import os
import numpy as np
import proposal as pp
from scipy.optimize import ridder
from importlib.resources import path

from taurunner.modules import units

def segment_body(body, granularity=0.5):
    descs = []
    for xi, xf in zip(body.layer_boundaries[1:], body.layer_boundaries[:-1]):
        if body.get_density(xf)/body.get_density(xi)>=granularity:
            descs.append((xi, xf, body.get_average_density(0.5*(xi+xf))))
        else:
            end   = 0
            start = xi
            while end<xf:
                s_density = body.get_density(start)
                func      = lambda x: body.get_density(x)-granularity*s_density
                end       = ridder(func, xi, xf)
                if end<xf:
                    I = quad(body.get_density, start, end, full_output=1)
                    avg_density  = I[0]/(end-start)
                    descs.append((start, end, avg_density))
                else:
                    I = quad(body.get_density, start, xf, full_output=1)
                    avg_density  = I[0]/(xf-start)
                    descs.append((start, xf, avg_density))
                start = end
    return descs

def make_propagator(body, xs_model='dipole', granularity=0.5):

    #define how many layers of constant density we need for the tau
    descs = segment_body(body, granularity)
    #make the sectors
    sec_defs = [make_sector(d/units.gr*units.cm**3, e*body.radius/units.cm, s*body.radius/units.cm, xs_model) for s, e, d in descs]
        
    with path('taurunner.resources.proposal_tables', 'tables.txt') as p:
        tables_path = str(p).split('tables.txt')[0]
    
    #define interpolator
    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = tables_path
    interpolation_def.path_to_tables_readonly = tables_path
    interpolation_def.nodes_cross_section = 200

    #define propagator -- takes a particle definition - sector - detector - interpolator
    prop = pp.Propagator(particle_def=pp.particle.TauMinusDef(),
                         sector_defs=sec_defs,
                         detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
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
    
    sec_def.cut_settings.ecut = 1e5*1e3
    sec_def.cut_settings.vcut = 0.1
    
    if(xs_model=='dipole'):
        sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.BlockDurandHa
    else:
        sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.AbramowiczLevinLevyMaor97
    
    return sec_def
