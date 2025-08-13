import numpy as np
import proposal as pp

from importlib.resources import path
from scipy.optimize import ridder
from scipy.integrate import quad
from typing import List, Tuple
from dataclasses import dataclass

from ..particle import Particle
from ..body import Body
from ..utils import units

PARTICLE_DEF_DICT = {
    11: pp.particle.EMinusDef,
    13: pp.particle.MuMinusDef,
    15: pp.particle.TauMinusDef,
    -11: pp.particle.EPlusDef,
    -13: pp.particle.MuPlusDef,
    -15: pp.particle.TauPlusDef,
}

@dataclass
class ConstantDensitySegment:
    ri: float
    rf: float
    density: float

def make_propagator(
    pdg_encoding: int,
    body: Body,
    xs_model
) -> pp.Propagator:
    """Make a PROPOSAL propagator

    params
    ______
    particle: Prometheus particle for which we want a PROPOSAL propagator
    simulation_specs: Dictionary specifying the configuration settings
    path_dict: Dictionary specifying any required path variables

    returns
    _______
                    make_propagatorprop: PROPOSAL propagator for input Particle
    """

    with path('taurunner.resources.proposal_tables', 'tables.txt') as p:
        tables_path = str(p).split('tables.txt')[0]
    #tables_path = "/n/home12/jlazar/TauRunner/taurunner/resources/proposal_tables/"
    pp.InterpolationSettings.tables_path = tables_path

    if pdg_encoding not in PARTICLE_DEF_DICT:
        raise ValueError(f"PDG code {pdg_encoding} is not a charged lepton")
    pdef = PARTICLE_DEF_DICT[pdg_encoding]()

    segments = segment_body(body)
    utilities = make_propagation_utilities(pdef, segments)
    geometries = make_geometries(segments, body)
    density_distrs = make_density_distributions(segments)
    prop = pp.Propagator(pdef, list(zip(geometries, utilities, density_distrs)))

    return prop

def make_geometries(segments: List[ConstantDensitySegment], body: Body) -> List[pp.Cartesian3D]:
    """Make list of proposal geometries from earth datafile

    params
    ______
    earth_file: data file where the parametrization of Earth is stored

    returns
    _______
    geometries: List of PROPOSAL spherical shells that make up the Earth
    """
    geometries = []
    for segment in segments:
        inner_radius = segment.ri * body.radius
        outer_radius = segment.rf * body.radius
        geometry = pp.geometry.Sphere(
            pp.Cartesian3D(0,0,1),
            outer_radius / units.cm,
            inner_radius / units.cm,
        )
        geometries.append(geometry)

    return geometries

def make_density_distributions(segments: List[ConstantDensitySegment]) -> List[pp.density_distribution.density_distribution]:
    """Make list of proposal homogeneous density distributions from
    Earth datafile

    params
    ______
    segments: A list of `ConstantDensitySegment`s

    returns
    _______
    density_distributions: Density distributions corresponding to the
        average density in each layer of the Earth model at linear order
    """
    density_distributions = []
    for segment in segments:
        density = pp.density_distribution.density_homogeneous(
            segment.density / (units.gr / units.cm**3)
        )
        density_distributions.append(density)
    return density_distributions

def segment_body(body: Body, granularity: float=0.5) -> List[Tuple[float, float, float]]:
    """
    Function for cutting a body into portions where the density does not grange to drastically.
    Since we use an approximation of constant density shells for propagation, this ensures that
    this approximation is not too crazy. Smaller granularity will make the approximation better
    at the expense of computation time.

    params
    ______
    body: TauRunner Body to be segmented
    granularity: fractional change in density at which to make a new layer

    returns
    _______
    segments: a list of tuples whose values are:
        0. the normalized radius at which a segment starts
        1. the normalized radius at which a segment ends
        2. the average density in this region
    """
    segments = []
    for xi, xf in zip(body.layer_boundaries[:-1], body.layer_boundaries[1:]):
        f_density = body.get_density(xf, right=True)
        gran = np.abs(f_density - body.get_density(xi, right=False)) / body.get_density(xi, right=False)
        if gran<granularity:  # the percent difference within a layer is small
            segments.append(
                ConstantDensitySegment(xi, xf, body.get_average_density(0.5 * (xi+xf)))
            )
            continue
        start = xi
        s_density = body.get_density(start, right=True)
        while np.abs(f_density-s_density)/s_density-granularity>0:
            func = lambda x: np.abs(body.get_density(x, right=True)-s_density)/s_density-granularity
            end = ridder(func, start, xf)
            wrap = lambda x: body.get_density(x, right=True)
            I = quad(wrap, start, end, full_output=1)
            avg_density  = I[0]/(end-start)
            segments.append(ConstantDensitySegment(start, end, avg_density))
            start = end
            s_density = body.get_density(start)
    return segments

def make_propagation_utilities(
    particle_def: pp.particle.ParticleDef,
    segments: List[ConstantDensitySegment],
) -> pp.PropagationUtility:
    """Make a list of PROPOSAL propagation utilities from an earth file
        for a particle given some simulation specifications

    params
    ______
    particle_def: PROPOSAL particle definition
    earth_file: data file where the parametrization of Earth is stored
    simulation_specs: dictionary specifying all the simulation specifications

    returns
    _______
    utilities: List of PROPOSAL PropagationUtility objects
    """
    cuts = pp.EnergyCutSettings(
        np.inf,
        1e-3,
        True
    )
    utilities = []
    components = [
        pp.component.Hydrogen(2),
        pp.component.Oxygen()
    ]
    for segment in segments:
        medium = pp.medium.Medium(
            f'Ice_{segment.density / (units.gr / units.cm**3)}',
            #1.0,
            75.0,
            -3.5017,
            0.09116,
            3.4773,
            0.2400,
            2.8004,
            0,
            segment.density / (units.gr / units.cm**3),
            components
        )
        collection = pp.PropagationUtilityCollection()
        cross = pp.crosssection.make_std_crosssection(
            particle_def,
            medium,
            cuts,
            True
        )
        create_tables = True
        collection.displacement = pp.make_displacement(cross, create_tables)
        collection.interaction = pp.make_interaction(cross, create_tables)
        collection.time = pp.make_time(cross, particle_def, create_tables)
        collection.decay = pp.make_decay(cross, particle_def, create_tables)

        utility = pp.PropagationUtility(collection=collection)
        utilities.append(utility)

    return utilities
