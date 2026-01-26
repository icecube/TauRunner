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

    if !PROPOSAL_AVAILABLE[]
        return nothing
    end

    # Create the appropriate propagator based on particle type
    pp = if particle_type == Tau
        @eval PROPOSAL.create_propagator_tauminus($(prop.config_path))
    elseif particle_type == AntiTau
        @eval PROPOSAL.create_propagator_tauplus($(prop.config_path))
    elseif particle_type == Muon
        @eval PROPOSAL.create_propagator_muminus($(prop.config_path))
    elseif particle_type == AntiMuon
        @eval PROPOSAL.create_propagator_muplus($(prop.config_path))
    elseif particle_type == Electron
        @eval PROPOSAL.create_propagator_eminus($(prop.config_path))
    elseif particle_type == Positron
        @eval PROPOSAL.create_propagator_eplus($(prop.config_path))
    else
        error("Unknown particle type for PROPOSAL: $particle_type")
    end

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

    # Use PROPOSAL if available, otherwise fall back to simplified model
    if PROPOSAL_AVAILABLE[]
        _propagate_with_proposal!(propagator, particle, track, remaining_dist_cm)
    else
        _simplified_propagation!(propagator, particle, track)
    end

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
    if isnothing(pp)
        _simplified_propagation!(propagator, particle, track)
        return nothing
    end

    # Get body radius in cm
    body_radius_cm = radius(propagator.body) / units.cm

    # Convert TauRunner position to PROPOSAL coordinates
    # TauRunner: x ∈ [0,1] along track, with x_to_r giving normalized radius
    # PROPOSAL: Cartesian coordinates in cm

    # Get current position in Cartesian coordinates
    # For a chord track, we need to compute the 3D position
    pos_cm, dir = track_to_proposal_coords(track, particle.position, body_radius_cm)

    # Energy in MeV for PROPOSAL
    energy_mev = particle.energy / units.MeV

    # Minimum energy (rest mass in MeV)
    min_energy_mev = particle.mass / units.MeV

    # Create PROPOSAL ParticleState
    particle_type_id = proposal_particle_type(particle.id)

    state = @eval PROPOSAL.ParticleState(
        $particle_type_id,
        $(pos_cm[1]), $(pos_cm[2]), $(pos_cm[3]),
        $(dir[1]), $(dir[2]), $(dir[3]),
        $energy_mev;
        time=0.0,
        propagated_distance=0.0
    )

    # Propagate
    secondaries = @eval PROPOSAL.propagate($pp, $state, $max_distance_cm, $min_energy_mev)

    # Extract final state
    final_state = @eval PROPOSAL.get_final_state($secondaries)
    final_energy_mev = @eval PROPOSAL.get_energy($final_state)
    propagated_dist_cm = @eval PROPOSAL.get_propagated_distance($final_state)

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
    if particle.energy <= particle.mass * 1.01  # Allow small tolerance
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
    # For a chord at angle θ starting at surface:
    # Position moves inward along the chord direction
    d = x_to_d(track, x)  # Normalized distance
    d_cm = d * body_radius_cm

    # Starting position (on surface, entry point)
    # For simplicity, assume entry at (0, 0, body_radius)
    # and direction pointing inward
    pos = (
        -dir[1] * d_cm,
        -dir[2] * d_cm,
        body_radius_cm - dir[3] * d_cm
    )

    return (pos, dir)
end

"""
Map TauRunner ParticleType to PROPOSAL particle type constant.
"""
function proposal_particle_type(p::ParticleType)
    # PROPOSAL particle type constants
    if p == Tau
        return @eval PROPOSAL.PARTICLE_TYPE_TAUMINUS
    elseif p == AntiTau
        return @eval PROPOSAL.PARTICLE_TYPE_TAUPLUS
    elseif p == Muon
        return @eval PROPOSAL.PARTICLE_TYPE_MUMINUS
    elseif p == AntiMuon
        return @eval PROPOSAL.PARTICLE_TYPE_MUPLUS
    elseif p == Electron
        return @eval PROPOSAL.PARTICLE_TYPE_EMINUS
    elseif p == Positron
        return @eval PROPOSAL.PARTICLE_TYPE_EPLUS
    else
        error("Unknown particle type: $p")
    end
end

# Import track coordinate functions
using ..TauRunner.Tracks: x_to_cartesian_direction

"""
Simplified propagation model (fallback when PROPOSAL.jl is not available).

⚠️  WARNING: This model is NOT suitable for physics analysis!

This is a rough approximation that:
1. Estimates energy loss based on path length and material (~2 MeV/(g/cm²))
2. Checks for decay based on mean decay length

Missing physics compared to PROPOSAL:
- Stochastic energy loss sampling (bremsstrahlung, pair production, photonuclear)
- Multiple scattering and angular deflection
- Proper Monte Carlo decay sampling
- Secondary particle production
"""
function _simplified_propagation!(
    propagator::SphericalBodyPropagator,
    particle,
    track
)
    # Warn user that simplified propagation is not physics-accurate
    warn_simplified_propagation()
    # Get remaining distance to exit (in natural units)
    remaining_x = 1.0 - particle.position
    remaining_dist = x_to_d(track, remaining_x) * radius(propagator.body)

    # Get average density along remaining path
    # (simplified: use density at current position)
    r = x_to_r(track, particle.position)
    density = get_density(propagator.body, r)

    # Rough energy loss estimate (for tau: ~2 MeV/(g/cm²) at high energies)
    # This is very approximate
    dE_dx = 2.0 * units.MeV * density / units.DENSITY_CONV

    column_depth = density * remaining_dist
    energy_loss = dE_dx * column_depth / (units.gr / units.cm^2)

    # Check if particle would decay before reaching exit
    # Decay length = γcτ = E/(mc²) * cτ
    if particle.lifetime < Inf && particle.mass > 0
        gamma = particle.energy / particle.mass
        decay_length = gamma * particle.lifetime

        if decay_length < remaining_dist
            # Particle decays before exit
            decay_frac = decay_length / remaining_dist
            particle.position = particle.position + decay_frac * remaining_x
            particle.decay_position = particle.position

            # Apply partial energy loss
            particle.energy = max(particle.energy - decay_frac * energy_loss, particle.mass)

            # Trigger decay
            decay!(particle)
            return nothing
        end
    end

    # Particle exits the body
    particle.position = 1.0
    particle.energy = max(particle.energy - energy_loss, particle.mass)

    return nothing
end
