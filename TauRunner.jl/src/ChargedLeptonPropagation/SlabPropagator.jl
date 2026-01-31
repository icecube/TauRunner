"""
Charged lepton propagation through slab bodies using PROPOSAL.jl.

Similar to SphericalBodyPropagator but for flat slab geometry.
"""

# Import Particles module for decay!
using ..TauRunner.Particles: Particle, decay!

"""
    SlabPropagator{B<:AbstractSlabBody} <: ChargedLeptonPropagator

Propagator for charged leptons through slab bodies.

# Fields
- `body::B`: The slab body
- `config_path::String`: Path to PROPOSAL JSON config file
- `propagators::Dict`: Cached PROPOSAL propagators by particle type
"""
mutable struct SlabPropagator{B<:AbstractSlabBody} <: ChargedLeptonPropagator
    body::B
    config_path::String
    propagators::Dict{ParticleType, Any}
end

"""
    SlabPropagator(body::AbstractSlabBody)

Create a propagator for a slab body.
"""
function SlabPropagator(body::B) where B<:AbstractSlabBody
    # Generate config from body
    config = generate_slab_config(body)
    config_path = create_temp_config(config)

    # Initialize propagator cache
    propagators = Dict{ParticleType, Any}()

    return SlabPropagator{B}(body, config_path, propagators)
end

"""
Get or create a PROPOSAL propagator for a given particle type.
"""
function get_proposal_propagator(prop::SlabPropagator, particle_type::ParticleType)
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
    propagate_charged_lepton!(propagator::SlabPropagator, particle, track)

Propagate a charged lepton through a slab body using PROPOSAL.jl.
"""
function propagate_charged_lepton!(
    propagator::SlabPropagator,
    particle,
    track::AbstractSlabTrack
)
    # Skip if losses are disabled or particle is electron
    abs_id = abs(Int(particle.id))
    if !particle.include_losses || abs_id == 11
        return nothing
    end

    # For neutrinos, nothing to do
    if abs_id in (12, 14, 16)
        return nothing
    end

    # Calculate remaining distance
    remaining_x = 1.0 - particle.position
    if remaining_x < 1e-10
        particle.position = 1.0
        return nothing
    end

    remaining_dist_cm = x_to_d(track, remaining_x) * propagator.body.length_natural / units.cm

    # For muons far from the exit (>100 km), they won't make it
    if abs_id == 13 && remaining_dist_cm / 1e5 > 100.0
        return nothing
    end

    _propagate_slab_with_proposal!(propagator, particle, track, remaining_dist_cm)

    return nothing
end

"""
Propagate through slab using PROPOSAL.jl.
"""
function _propagate_slab_with_proposal!(
    propagator::SlabPropagator,
    particle,
    track,
    max_distance_cm::Float64
)
    pp = get_proposal_propagator(propagator, particle.id)

    # Slab length in cm
    slab_length_cm = propagator.body.length_natural / units.cm

    # Current position in slab (in cm)
    # Clamp to stay strictly inside the geometry — PROPOSAL segfaults if the
    # particle position lands exactly on or beyond the slab boundary due to
    # floating-point arithmetic ("No sector defined at particle position").
    current_z_cm = clamp(particle.position * slab_length_cm, 0.0, slab_length_cm * (1.0 - 1e-12))

    # Direction (along z-axis for normal incidence, or tilted for angled tracks)
    dir = x_to_cartesian_direction(track, particle.position)

    # Energy in MeV
    energy_mev = particle.energy / units.MeV
    # Minimum energy: use 0 to match Python TauRunner (PROPOSAL default)
    min_energy_mev = 0.0

    # Create PROPOSAL ParticleState
    particle_type_id = proposal_particle_type(particle.id)

    state = Base.invokelatest(
        PROPOSAL_FN[:ParticleState],
        particle_type_id,
        0.0, 0.0, current_z_cm,
        dir[1], dir[2], dir[3],
        energy_mev;
        time=0.0,
        propagated_distance=0.0
    )

    # Propagate
    secondaries = suppress_proposal_warnings() do
        Base.invokelatest(PROPOSAL_FN[:propagate], pp, state, max_distance_cm, min_energy_mev)
    end

    # Extract final state
    final_state = Base.invokelatest(PROPOSAL_FN[:get_final_state], secondaries)
    final_energy_mev = Base.invokelatest(PROPOSAL_FN[:get_energy], final_state)
    propagated_dist_cm = Base.invokelatest(PROPOSAL_FN[:get_propagated_distance], final_state)

    # Update particle energy
    particle.energy = final_energy_mev * units.MeV

    # Update position
    dist_traveled_normalized = (propagated_dist_cm * units.cm) / propagator.body.length_natural
    particle.position = clamp(particle.position + dist_traveled_normalized, 0.0, 1.0)

    # Check for decay
    if particle.energy <= particle.mass * (1.0 + 1e-6)
        particle.decay_position = particle.position
        decay!(particle)
    end

    return nothing
end

# Import track functions
using ..TauRunner.Tracks: x_to_cartesian_direction
