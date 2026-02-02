"""
Charged lepton propagation through slab bodies using PROPOSAL.jl.

Matches the Python TauRunner approach: each layer gets its own propagator
with a large enclosing sphere geometry, and the particle always starts at
the origin. This avoids PROPOSAL segfaults from boundary issues.
"""

# Import Particles module for decay!
using ..TauRunner.Particles: Particle, decay!

"""
    SlabPropagator{B<:AbstractSlabBody} <: ChargedLeptonPropagator

Propagator for charged leptons through slab bodies.

# Fields
- `body::B`: The slab body
- `propagators::Dict`: Cached PROPOSAL propagators by (particle_type, medium_name)
"""
mutable struct SlabPropagator{B<:AbstractSlabBody} <: ChargedLeptonPropagator
    body::B
    propagators::Dict{Tuple{ParticleType, String}, Any}
end

"""
    SlabPropagator(body::AbstractSlabBody)

Create a propagator for a slab body.
"""
function SlabPropagator(body::B) where B<:AbstractSlabBody
    propagators = Dict{Tuple{ParticleType, String}, Any}()
    return SlabPropagator{B}(body, propagators)
end

"""
Get or create a PROPOSAL propagator for a given particle type and medium.
Each (particle_type, medium) pair gets its own propagator with a large
enclosing sphere, matching the Python TauRunner approach.
"""
function get_proposal_propagator(prop::SlabPropagator, particle_type::ParticleType, medium_name::String, density_gcm3::Float64)
    key = (particle_type, medium_name)
    if haskey(prop.propagators, key)
        return prop.propagators[key]
    end

    # Generate config for this layer
    config = generate_slab_layer_config(density_gcm3, medium_name)
    config_path = create_temp_config(config)

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
    pp = Base.invokelatest(PROPOSAL_FN[creator_key], config_path)

    prop.propagators[key] = pp
    return pp
end

"""
    propagate_charged_lepton!(propagator::SlabPropagator, particle, track)

Propagate a charged lepton through a slab body using PROPOSAL.jl.
Iterates over layers like the Python implementation.
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

    # For muons far from the exit (>100 km), they won't make it
    remaining_dist_cm = x_to_d(track, remaining_x) * propagator.body.length_natural / units.cm
    if abs_id == 13 && remaining_dist_cm / 1e5 > 100.0
        return nothing
    end

    _propagate_slab_with_proposal!(propagator, particle, track)

    return nothing
end

"""
Propagate through slab layer-by-layer using PROPOSAL.jl.

Matches the Python TauRunner approach: for each layer that the particle
still needs to traverse, create a propagator with a large enclosing sphere
and start the particle at the origin. This completely avoids PROPOSAL
segfaults from particles landing on sector boundaries.
"""
function _propagate_slab_with_proposal!(
    propagator::SlabPropagator,
    particle,
    track
)
    body = propagator.body
    boundaries = layer_boundaries(body)
    slab_length_cm = body.length_natural / units.cm

    particle_type_id = proposal_particle_type(particle.id)
    dir = x_to_cartesian_direction(track, particle.position)
    min_energy_mev = 0.0

    prv = particle.position

    for i in 1:(length(boundaries) - 1)
        layer_end = boundaries[i + 1]

        # Skip layers the particle has already passed
        if layer_end <= particle.position
            continue
        end

        # Get layer properties
        x_mid = (boundaries[i] + boundaries[i + 1]) / 2
        density = get_density(body, x_mid)
        density_gcm3 = density / units.DENSITY_CONV
        medium_name = density_to_medium(density)

        # Get or create propagator for this (particle_type, medium) pair
        pp = get_proposal_propagator(propagator, particle.id, medium_name, density_gcm3)

        # Distance to propagate in this layer (in cm)
        prop_length_cm = (layer_end - prv) * slab_length_cm

        # Energy in MeV
        energy_mev = particle.energy / units.MeV

        # Create PROPOSAL ParticleState at the origin (like Python)
        state = Base.invokelatest(
            PROPOSAL_FN[:ParticleState],
            particle_type_id,
            0.0, 0.0, 0.0,
            dir[1], dir[2], dir[3],
            energy_mev;
            time=0.0,
            propagated_distance=0.0
        )

        # Propagate for the layer thickness
        secondaries = suppress_proposal_warnings() do
            Base.invokelatest(PROPOSAL_FN[:propagate], pp, state, prop_length_cm, min_energy_mev)
        end

        # Extract final state
        final_state = Base.invokelatest(PROPOSAL_FN[:get_final_state], secondaries)
        final_energy_mev = Base.invokelatest(PROPOSAL_FN[:get_energy], final_state)
        propagated_dist_cm = Base.invokelatest(PROPOSAL_FN[:get_propagated_distance], final_state)

        # Update particle energy
        particle.energy = final_energy_mev * units.MeV

        # Update position based on actual distance propagated
        dist_traveled_normalized = (propagated_dist_cm * units.cm) / body.length_natural
        particle.position = clamp(prv + dist_traveled_normalized, 0.0, 1.0)

        # Check for decay
        if particle.energy <= particle.mass * (1.0 + 1e-6)
            particle.decay_position = particle.position
            decay!(particle)
            return nothing
        end

        prv = layer_end
    end

    return nothing
end

# Import track functions
using ..TauRunner.Tracks: x_to_cartesian_direction
