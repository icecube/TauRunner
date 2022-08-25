import numpy as np
import proposal as pp
from scipy.optimize import ridder
import sys
# In Python 3.6 and before importlib.resources is importlib_resources
if sys.version_info.major==3 and sys.version_info.minor<=6:
    from importlib_resources import path
else:
    from importlib.resources import path

from scipy.integrate import quad
from taurunner.utils import units

ID_2_name = {
    16: 'TauMinusDef',
    -16: 'TauPlusDef',
    14: 'MuMinusDef',
    -14: 'MuPlusDef',
    12: 'EMinusDef',
    -12:'EPlusDef'
}

def segment_body(body, granularity=0.5):
    descs = []
    for xi, xf in zip(body.layer_boundaries[:-1], body.layer_boundaries[1:]):
        f_density = body.get_density(xf, right=True)
        gran = np.abs(f_density - body.get_density(xi, right=False)) \
            / body.get_density(xi, right=False)
        if gran<granularity:  # the percent difference within a layer is small
            descs.append((xi, xf, body.get_average_density(0.5*(xi+xf))))
        else:
            start = xi
            s_density = body.get_density(start, right=True)
            while np.abs(f_density-s_density)/s_density-granularity>0:
                func = lambda x: np.abs(body.get_density(x, right=True)-s_density)/s_density-granularity
                end = ridder(func, start, xf)
                wrap = lambda x: body.get_density(x, right=True)
                I = quad(wrap, start, end, full_output=1)
                avg_density  = I[0]/(end-start)
                descs.append((start, end, avg_density))
                start = end
                s_density = body.get_density(start)
    return descs

def make_propagator(ID, body, xs_model='dipole', granularity=0.5):
    
    if(ID in [12, -12]):
        return None

    #with path('taurunner.resources.proposal_tables', 'tables.txt') as p:
    #    pp.InterpolationSettings.tables_path = str(p).split('tables.txt')[0]

    # Make the photonuclear modeldescs = segment_body(body, granularity)
    shadow_effect = pp.parametrization.photonuclear.ShadowDuttaRenoSarcevicSeckel()
    if xs_model=="dipole":
       parametrization = (
            pp.parametrization.photonuclear.BlockDurandHa(shadow_effect)
        )
    else:
        parametrization = (
            pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97(shadow_effect)
        )
    # Make the particle
    particle_def = getattr(pp.particle, ID_2_name[ID])()
    # Configure energy loss settings
    ecut = 200
    vcut = 1e-3
    do_continuous_randomization = True
    interpolate = True
    cuts = pp.EnergyCutSettings(ecut, vcut, do_continuous_randomization)
    #define layers of constant density we need for the tau
    descs = segment_body(body, granularity)
    # Hacky way to add air padding...
    descs.append((1, 1.2, 0.001))
    # make the PROPOSAL propagation utilities
    utilities = [make_utility(
        particle_def,
        interpolate,
        cuts,
        parametrization,
        density/units.gr*units.cm**3
    ) for _, _, density in descs]

    # Make the PROPOSAL geometries
    geometries = [
        pp.geometry.Sphere(
            pp.Cartesian3D(),
            end*body.radius/units.cm,
            start*body.radius/units.cm
        )
        for start, end, _ in descs
    ]
    # Make the PROPOSAL density distributions
    density_distrs = [
        pp.density_distribution.density_homogeneous(
            density/units.gr*units.cm**3
        ) for _, _, density in descs
    ]

    prop = pp.Propagator(
        particle_def,
        list(zip(geometries, utilities, density_distrs))
    )

    return prop

def make_utility(
    particle_def,
    interpolate,
    cuts,
    parametrization,
    density
):
    collection = pp.PropagationUtilityCollection()
    components = [
        pp.component.Hydrogen(2),
        pp.component.Oxygen()
    ]
    target = pp.medium.Medium(
        f'Ice_{density}',
        75.0,
        -3.5017, 
        0.09116,
        3.4773,
        0.2400, 
        2.8004,
        0,
        density,
        components
    ) 
    cross = pp.crosssection.make_std_crosssection(
        particle_def,
        target,
        cuts,
        interpolate
    )
    cross[3] = pp.crosssection.make_crosssection(
        parametrization,
        particle_def,
        target,
        cuts,
        interpolate
    )
    # TODO look into continuous randomization
    create_tables = True
    collection.displacement = pp.make_displacement(cross, create_tables)
    collection.interaction = pp.make_interaction(cross, create_tables)
    collection.time = pp.make_time(cross, particle_def, create_tables)
    collection.decay = pp.make_decay(cross, particle_def, create_tables)
    
    utility = pp.PropagationUtility(collection=collection)

    return utility

def make_sector(density, start, end, xs_model):
    #Define a sector

    #sec_def = pp.SectorDefinition()
    components = [
        pp.component.Hydrogen(2),
        pp.component.Oxygen()
    ]
    medium = pp.medium.Medium(
        f'Ice_{density}',
        75.0,
        -3.5017, 
        0.09116,
        3.4773,
        0.2400, 
        2.8004,
        0,
        density,
        components
    ) 
    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), end, start)
    sec_def.particle_location = pp.ParticleLocation.inside_detector        
    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere
    sec_def.crosssection_defs.brems_def.lpm_effect = True
    sec_def.crosssection_defs.epair_def.lpm_effect = True
    sec_def.cut_settings.ecut = -1.0
    sec_def.cut_settings.vcut = 1e-3
    sec_def.do_continuous_randomization = True

    if(xs_model=='dipole'):
        sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.BlockDurandHa
    else:
        sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.AbramowiczLevinLevyMaor97

    return sec_def
