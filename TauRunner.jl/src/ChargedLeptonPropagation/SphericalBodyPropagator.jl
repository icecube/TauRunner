"""
Charged lepton propagation through spherical bodies using PROPOSAL.jl.

This file provides the interface to PROPOSAL.jl for propagating charged leptons
(tau, muon) through spherical bodies like Earth or Sun.
"""

# Import Particles module for decay!
using ..TauRunner.Particles: Particle, decay!

"""
    SphericalBodyPropagator{B<:AbstractSphericalBody} <: ChargedLeptonPropagator

Propagator for charged leptons through spherical bodies.

# Fields
- `body::B`: The spherical body
- `config_path::String`: Path to PROPOSAL JSON config file
- `propagators::Dict`: Cached PROPOSAL propagators by particle type
"""
mutable struct SphericalBodyPropagator{B<:AbstractSphericalBody} <: ChargedLeptonPropagator
    body::B
    config_path::String
    propagators::Dict{ParticleType, Any}
end

"""
    SphericalBodyPropagator(body::AbstractSphericalBody)

Create a propagator for a spherical body.

Generates a PROPOSAL configuration based on the body's layer structure
and initializes PROPOSAL propagators for each charged lepton type.
"""
function SphericalBodyPropagator(body::B) where B<:AbstractSphericalBody
    # Generate config from body
    config = generate_sphere_config(body)
    config_path = create_temp_config(config)

    # Initialize propagator cache (will be populated on first use)
    propagators = Dict{ParticleType, Any}()

    return SphericalBodyPropagator{B}(body, config_path, propagators)
end

"""
Get or create a PROPOSAL propagator for a given particle type.
"""
function get_proposal_propagator(prop::SphericalBodyPropagator, particle_type::ParticleType)
    if haskey(prop.propagators, particle_type)
        return prop.propagators[particle_type]
    end

    # Create the appropriate propagator based on particle type
    creator_key = if particle_type == Tau
        :create_propagator_tauminus
    elseif particle_type == AntiTau
        :create_propagator_tauplus
    elseif particle_type == Muon
        :create_propagator_muminus
    elseif particle_type == AntiMuon
        :create_propagator_muplus
    elseif particle_type == Electron
        :create_propagator_eminus
    elseif particle_type == Positron
        :create_propagator_eplus
    else
        error("Unknown particle type for PROPOSAL: $particle_type")
    end
    pp = Base.invokelatest(PROPOSAL_FN[creator_key], prop.config_path)

    prop.propagators[particle_type] = pp
    return pp
end

"""
    propagate_charged_lepton!(propagator::SphericalBodyPropagator, particle, track)

Propagate a charged lepton through a spherical body using PROPOSAL.jl.

# Arguments
- `propagator`: SphericalBodyPropagator instance
- `particle`: TauRunner Particle to propagate (modified in place)
- `track`: Track geometry (Chord, Radial, etc.)

Updates the particle's position, energy, and handles decay if applicable.
"""
function propagate_charged_lepton!(
    propagator::SphericalBodyPropagator,
    particle,
    track::AbstractSphericalTrack
)
    # Skip if losses are disabled or particle is electron (immediately absorbed)
    abs_id = abs(Int(particle.id))
    if !particle.include_losses || abs_id == 11
        return nothing
    end

    # For neutrinos, nothing to do
    if abs_id in (12, 14, 16)
        return nothing
    end

    # Calculate remaining distance to exit
    remaining_x = 1.0 - particle.position
    if remaining_x < 1e-10
        particle.position = 1.0
        return nothing
    end

    # Distance in natural units, then convert to cm for PROPOSAL
    remaining_dist_natural = x_to_d(track, remaining_x) * radius(propagator.body)
    remaining_dist_cm = remaining_dist_natural / units.cm

    # For muons far from the exit (>100 km), they won't make it
    if abs_id == 13 && remaining_dist_cm / 1e5 > 100.0  # 100 km in cm
        return nothing
    end

    _propagate_with_proposal!(propagator, particle, track, remaining_dist_cm)

    return nothing
end

"""
Propagate using PROPOSAL.jl.
"""
function _propagate_with_proposal!(
    propagator::SphericalBodyPropagator,
    particle,
    track,
    max_distance_cm::Float64
)
    pp = get_proposal_propagator(propagator, particle.id)

    # Get body radius in cm
    body_radius_cm = radius(propagator.body) / units.cm

    # Convert TauRunner position to PROPOSAL coordinates
    # TauRunner: x ∈ [0,1] along track, with x_to_r giving normalized radius
    # PROPOSAL: Cartesian coordinates in cm

    # Get current position in Cartesian coordinates
    # For a chord track, we need to compute the 3D position
    # Clamp position to stay strictly inside the geometry — PROPOSAL errors if the
    # particle position lands exactly on or beyond a boundary due to floating-point
    # arithmetic ("No sector defined at particle position").
    clamped_x = clamp(particle.position, 1e-12, 1.0 - 1e-12)
    pos_cm, dir = track_to_proposal_coords(track, clamped_x, body_radius_cm)

    # Energy in MeV for PROPOSAL
    energy_mev = particle.energy / units.MeV

    # Minimum energy: use 0 to match Python TauRunner (PROPOSAL default)
    min_energy_mev = 0.0

    # Create PROPOSAL ParticleState
    particle_type_id = proposal_particle_type(particle.id)

    state = Base.invokelatest(
        PROPOSAL_FN[:ParticleState],
        particle_type_id,
        pos_cm[1], pos_cm[2], pos_cm[3],
        dir[1], dir[2], dir[3],
        energy_mev;
        time=0.0,
        propagated_distance=0.0
    )

    # Propagate
    # Round distance to integer cm to avoid PROPOSAL bisection failures at geometry boundaries
    # (matches Python TauRunner convention)
    max_distance_cm_int = floor(max_distance_cm)

    secondaries = suppress_proposal_warnings() do
        Base.invokelatest(PROPOSAL_FN[:propagate], pp, state, max_distance_cm_int, min_energy_mev)
    end

    # Extract final state
    final_state = Base.invokelatest(PROPOSAL_FN[:get_final_state], secondaries)
    final_energy_mev = Base.invokelatest(PROPOSAL_FN[:get_energy], final_state)
    propagated_dist_cm = Base.invokelatest(PROPOSAL_FN[:get_propagated_distance], final_state)

    # Update particle energy
    particle.energy = final_energy_mev * units.MeV

    # Update particle position based on propagated distance
    dist_traveled_natural = propagated_dist_cm * units.cm
    dist_traveled_normalized = dist_traveled_natural / radius(propagator.body)

    # Convert distance traveled to new track parameter
    current_d = x_to_d(track, particle.position)
    new_d = current_d + dist_traveled_normalized
    new_x = d_to_x(track, new_d)

    # Clamp to [0, 1]
    particle.position = clamp(new_x, 0.0, 1.0)

    # Check if particle decayed (energy dropped to rest mass or below)
    # For tau/muon, if energy is at rest mass, it decayed
    if particle.energy <= particle.mass * (1.0 + 1e-6)  # Rest mass = decayed
        particle.decay_position = particle.position
        decay!(particle)
    end

    return nothing
end

"""
Convert track position to PROPOSAL Cartesian coordinates.

Returns (position_cm, direction_unit_vector).
"""
function track_to_proposal_coords(track::AbstractSphericalTrack, x::Real, body_radius_cm::Real)
    # Get direction from track
    dir = x_to_cartesian_direction(track, x)

    # Get radius at this position
    r = x_to_r(track, x)
    r_cm = r * body_radius_cm

    # Compute position from track geometry
    # Use remaining distance (1-x) to place the particle correctly:
    # At x=0 (entry), the particle is at the entry point on the surface
    # At x=1 (exit), the particle is at (0, 0, body_radius)
    # This matches Python TauRunner's x_to_pp_pos convention
    remaining_d = x_to_d(track, 1.0) - x_to_d(track, x)  # Remaining normalized distance to exit
    remaining_d_cm = remaining_d * body_radius_cm

    pos = (
        -dir[1] * remaining_d_cm,
        -dir[2] * remaining_d_cm,
        (1.0 - track.depth) * body_radius_cm - dir[3] * remaining_d_cm
    )

    return (pos, dir)
end

"""
Map TauRunner ParticleType to PROPOSAL particle type constant.
"""
function proposal_particle_type(p::ParticleType)
    if p == Tau
        return PROPOSAL_FN[:PARTICLE_TYPE_TAUMINUS]
    elseif p == AntiTau
        return PROPOSAL_FN[:PARTICLE_TYPE_TAUPLUS]
    elseif p == Muon
        return PROPOSAL_FN[:PARTICLE_TYPE_MUMINUS]
    elseif p == AntiMuon
        return PROPOSAL_FN[:PARTICLE_TYPE_MUPLUS]
    elseif p == Electron
        return PROPOSAL_FN[:PARTICLE_TYPE_EMINUS]
    elseif p == Positron
        return PROPOSAL_FN[:PARTICLE_TYPE_EPLUS]
    else
        error("Unknown particle type: $p")
    end
end

# Import track coordinate functions
using ..TauRunner.Tracks: x_to_cartesian_direction
