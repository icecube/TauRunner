"""
Particle state management for Monte Carlo simulation.
"""

# Isoscalar mass (average of proton and neutron masses)
const ISOSCALAR_MASS = PhysicsConstants.ISOSCALAR_MASS_GEV * units.GeV

# Minimum energy allowed in cross-section splines
const EMIN = 1e9  # eV

# Energy fractions for sampling differential cross-sections
const NEUTRINO_DIFF_ENERGY_FRACTIONS = collect(range(0.0, 1.0, length=101))

"""
    SecondaryParticle{T}

A secondary particle produced in a decay.

# Fields
- `id::ParticleType`: PDG particle identifier
- `energy::T`: Energy in eV
- `position::T`: Track position where produced (0 to 1)
"""
struct SecondaryParticle{T<:Real}
    id::ParticleType
    energy::T
    position::T
end

"""
    Particle{T, XS}

A particle being propagated through matter.

# Fields
- `id::ParticleType`: PDG particle identifier
- `energy::T`: Current energy in eV
- `position::T`: Current position along track (0 to 1)
- `initial_energy::T`: Starting energy (for bookkeeping)
- `survived::Bool`: Whether particle has survived (not absorbed/stopped)
- `xs::XS`: Cross-sections object
- `include_secondaries::Bool`: Whether to track secondary particles
- `include_losses::Bool`: Whether to include energy losses
- `n_cc::Int`: Number of CC interactions
- `n_nc::Int`: Number of NC interactions
- `n_decay::Int`: Number of decays
- `decay_position::T`: Position of last decay
- `basket::Vector{SecondaryParticle{T}}`: Secondary particles from decays
"""
mutable struct Particle{T<:Real, XS<:AbstractCrossSections}
    id::ParticleType
    energy::T
    position::T
    initial_energy::T
    survived::Bool
    xs::XS
    include_secondaries::Bool
    include_losses::Bool
    n_cc::Int
    n_nc::Int
    n_decay::Int
    decay_position::T
    basket::Vector{SecondaryParticle{T}}
    # Particle properties (cached)
    mass::T
    lifetime::T
end

"""
    Particle(id, energy, position, xs; secondaries=true, losses=true)

Create a new particle.

# Arguments
- `id`: Particle type (ParticleType enum or integer PDG code)
- `energy`: Initial energy in eV
- `position`: Initial track position (0 to 1)
- `xs`: CrossSections object
- `secondaries`: Whether to track secondaries from decays (default: true)
- `losses`: Whether to include energy losses (default: true)
"""
function Particle(
    id::Union{ParticleType, Integer},
    energy::Real,
    position::Real,
    xs::XS;
    secondaries::Bool=true,
    losses::Bool=true
) where XS<:AbstractCrossSections
    T = Float64

    # Convert integer to ParticleType if needed
    particle_id = id isa Integer ? ParticleType(id) : id

    # Get particle properties
    mass, lifetime = get_particle_properties(particle_id)

    return Particle{T, XS}(
        particle_id,
        T(energy),
        T(position),
        T(energy),       # initial_energy
        true,            # survived
        xs,
        secondaries,
        losses,
        0,               # n_cc
        0,               # n_nc
        0,               # n_decay
        T(0),            # decay_position
        SecondaryParticle{T}[],
        T(mass),
        T(lifetime)
    )
end

"""
Get mass and lifetime for a particle type.
"""
function get_particle_properties(id::ParticleType)
    abs_id = abs(Int(id))

    if abs_id in (12, 14, 16)  # Neutrinos
        return (0.0, Inf)
    elseif abs_id == 15  # Tau
        return (PhysicsConstants.TAU_MASS_GEV * units.GeV,
                PhysicsConstants.TAU_LIFETIME_SEC * units.sec)
    elseif abs_id == 13  # Muon
        return (PhysicsConstants.MUON_MASS_GEV * units.GeV,
                PhysicsConstants.MUON_LIFETIME_SEC * units.sec)
    elseif abs_id == 11  # Electron
        return (PhysicsConstants.ELECTRON_MASS_GEV * units.GeV, Inf)
    else
        error("Unknown particle type: $id")
    end
end

"""
    nutype(particle::Particle)

Get the neutrino type string for cross-section lookups.
"""
function nutype(particle::Particle)
    return Int(particle.id) < 0 ? "nubar" : "nu"
end

"""
    get_proposed_depth_step(particle::Particle; rng=Random.default_rng())

Sample the free-streaming column depth before next interaction.

Returns the column depth to next interaction in natural units.
"""
function get_proposed_depth_step(particle::Particle{T}; rng::AbstractRNG=Random.default_rng()) where T
    if !is_neutrino(particle.id)
        error("get_proposed_depth_step is only valid for neutrinos")
    end

    p = rand(rng)

    cc_depth = get_interaction_depth(particle, :CC)
    nc_depth = get_interaction_depth(particle, :NC)

    # Combined interaction rate
    total_rate = 1/cc_depth + 1/nc_depth

    # Sample from exponential distribution
    return T(-log(p) / total_rate)
end

"""
    get_interaction_depth(particle::Particle, interaction; proton_fraction=0.5)

Calculate the mean column depth to interaction.

# Arguments
- `particle`: The particle
- `interaction`: :CC or :NC
- `proton_fraction`: Fraction of target that is protons

# Returns
Mean interaction depth in natural units
"""
function get_interaction_depth(
    particle::Particle{T},
    interaction::Symbol;
    proton_fraction::Real=0.5
) where T
    if !is_neutrino(particle.id)
        error("get_interaction_depth is only valid for neutrinos")
    end

    sigma = total_cross_section(
        particle.xs,
        particle.energy,
        nutype(particle),
        String(interaction),
        proton_fraction=proton_fraction
    )

    return T(ISOSCALAR_MASS / sigma)
end

"""
    get_total_interaction_depth(particle::Particle)

Get the combined interaction depth for all interactions.
"""
function get_total_interaction_depth(particle::Particle{T}) where T
    cc_depth = get_interaction_depth(particle, :CC)
    nc_depth = get_interaction_depth(particle, :NC)
    return T(1 / (1/cc_depth + 1/nc_depth))
end

"""
    interact!(particle::Particle, interaction; rng=Random.default_rng(), proton_fraction=0.5)

Perform a neutrino interaction, updating particle state.

# Arguments
- `particle`: The particle (modified in place)
- `interaction`: :CC or :NC
- `rng`: Random number generator
- `proton_fraction`: Fraction of target that is protons
"""
function interact!(
    particle::Particle{T},
    interaction::Symbol;
    rng::AbstractRNG=Random.default_rng(),
    proton_fraction::Real=0.5
) where T
    if !is_neutrino(particle.id)
        error("interact! is only valid for neutrinos, got $(particle.id)")
    end

    # Sample energy fraction from differential cross-section
    weights = differential_cross_section(
        particle.xs,
        particle.energy,
        NEUTRINO_DIFF_ENERGY_FRACTIONS,
        nutype(particle),
        String(interaction),
        proton_fraction=proton_fraction
    )

    # Normalize weights in-place to avoid allocation
    w_sum = sum(weights)
    weights ./= w_sum

    # Sample energy fraction
    z_choice = sample_weighted(rng, NEUTRINO_DIFF_ENERGY_FRACTIONS, weights)

    # Update energy (with minimum energy cutoff)
    particle.energy = z_choice * (particle.energy - EMIN) + EMIN

    if interaction == :CC
        # Charged current: neutrino -> charged lepton
        particle.n_cc += 1
        particle.id = charged_partner(particle.id)
        particle.mass, particle.lifetime = get_particle_properties(particle.id)

        # Electrons are immediately absorbed
        if abs(Int(particle.id)) == 11
            particle.survived = false
        end
    elseif interaction == :NC
        # Neutral current: neutrino stays neutrino
        particle.n_nc += 1
    else
        error("Unknown interaction type: $interaction")
    end

    return nothing
end

"""
Sample from a weighted distribution.

Uses a running cumulative sum to avoid allocating a temporary array.
"""
function sample_weighted(rng::AbstractRNG, values::AbstractVector, weights::AbstractVector)
    r = rand(rng)
    cumulative = zero(eltype(weights))
    @inbounds for i in eachindex(weights)
        cumulative += weights[i]
        if cumulative >= r
            return values[i]
        end
    end
    return @inbounds values[end]
end

"""
    reset!(particle::Particle, energy, position)

Reset particle state for reuse.
"""
function reset!(particle::Particle{T}, energy::Real, position::Real) where T
    particle.energy = T(energy)
    particle.position = T(position)
    particle.initial_energy = T(energy)
    particle.survived = true
    particle.n_cc = 0
    particle.n_nc = 0
    particle.n_decay = 0
    particle.decay_position = T(0)
    empty!(particle.basket)
    return nothing
end

function Base.show(io::IO, p::Particle)
    E_GeV = p.energy / units.GeV
    print(io, "Particle($(p.id), E=$(round(E_GeV, sigdigits=3)) GeV, x=$(round(p.position, digits=3)))")
end

# Cross-section functions are imported by the parent module (Particles.jl)
