import os
import numpy as np
import proposal as pp


def make_propagator():
    #Make proposal object
    #define geometry
    density = 1.0
    sec_def = make_sector(density)
       
    #define interpolator
    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "/home/isafa/.local/share/PROPOSAL/tables"
    interpolation_def.path_to_tables_readonly = "/home/isafa/.local/share/PROPOSAL/tables"
    interpolation_def.nodes_cross_section = 200
        
    #define propagator -- takes a particle definition - sector - detector - interpolator
    prop = pp.Propagator(particle_def=pp.particle.TauMinusDef(),
                         sector_defs=[sec_def],
                         detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
                         interpolation_def=interpolation_def)

    return prop


def make_sector(density):
    sec_def = pp.SectorDefinition()
    sec_def.medium = pp.medium.Ice(1.0)
    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), 1e20, 0)
    sec_def.particle_location = pp.ParticleLocation.inside_detector        
    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere
    sec_def.crosssection_defs.brems_def.lpm_effect = True
    sec_def.crosssection_defs.epair_def.lpm_effect = True
    
    sec_def.cut_settings.ecut = 1e5*1e3
    sec_def.cut_settings.vcut = 0.1

    sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.BlockDurandHa
    
    return sec_def
