"""
Casino.jl - Core Monte Carlo propagation algorithm for TauRunner.jl

This is the main propagation engine that simulates particle transport through matter.
Named "Casino" because Monte Carlo methods are central to the simulation.
"""

using Random

# Import from submodules
using .Bodies: get_density, radius, length
using .Tracks: total_column_depth, X_to_x, x_to_X, x_to_d
using .CrossSectionsModule: CrossSections, total_cross_section
using .Particles: Particle, SecondaryParticle, get_proposed_depth_step, get_interaction_depth
using .Particles: get_total_interaction_depth, interact!, decay!, reset!
using .ChargedLeptonPropagationModule: ChargedLeptonPropagator, propagate_charged_lepton!
using .ChargedLeptonPropagationModule: SphericalBodyPropagator, SlabPropagator

"""
    propagate!(particle, track, body, clp; condition=nothing, rng=Random.default_rng())

Simulate particle propagation along a track through a body.

This is the main Monte Carlo propagation function. It handles:
- Neutrino interactions (CC and NC)
- Charged lepton propagation with energy losses
- Particle decays
- Secondary particle tracking

# Arguments
- `particle`: The particle to propagate (modified in place)
- `track`: Track geometry (Chord, Radial, SlabTrack, etc.)
- `body`: Body through which to propagate (Earth, Sun, slab, etc.)
- `clp`: ChargedLeptonPropagator for handling tau/muon propagation
- `condition`: Optional stopping condition function `(particle, body, track) -> Bool`
- `rng`: Random number generator

# Returns
The particle after propagation (same object, modified in place)

# Algorithm
1. For neutrinos:
   - Sample interaction depth from exponential distribution
   - Advance position until interaction or exit
   - Sample interaction type (CC or NC) based on cross-sections
   - Apply interaction (energy loss, possible flavor change)
   - If CC: create charged lepton and propagate
2. For charged leptons:
   - Propagate with energy losses using PROPOSAL
   - Handle decay → back to tau neutrino

# Example
```julia
earth = construct_earth()
xs = CrossSections(CSMS)
track = Chord(π/4)  # 45 degree angle
clp = SphericalBodyPropagator(earth)

particle = Particle(NuTau, 1e18, 0.0, xs)
result = propagate!(particle, track, earth, clp)
```
"""
function propagate!(
    particle::Particle,
    track::AbstractTrack,
    body::AbstractBody,
    clp::ChargedLeptonPropagator;
    condition=nothing,
    rng::AbstractRNG=Random.default_rng()
)
    # Get total column depth for this track/body combination
    total_depth = total_column_depth(track, body)

    accumulated_depth = 0.0

    # Main propagation loop
    # Inline default stop check to avoid closure overhead on every iteration
    while true
        # Default stop: exited or absorbed
        (particle.position >= 1.0 || !particle.survived) && break
        # User-provided stop condition
        (!isnothing(condition) && condition(particle, body, track)) && break
        if is_neutrino(particle.id)
            # Neutrino propagation: sample interaction
            depth_step = get_proposed_depth_step(particle; rng=rng)
            accumulated_depth += depth_step

            # Check if we've passed the total column depth
            if accumulated_depth >= total_depth
                particle.position = 1.0
                return particle
            end

            # Update position based on accumulated column depth
            particle.position = X_to_x(track, body, accumulated_depth)

            # Determine interaction type
            p_rand = rand(rng)
            cc_depth = get_interaction_depth(particle, :CC)
            total_int_depth = get_total_interaction_depth(particle)
            p_cc = total_int_depth / cc_depth

            interaction = p_rand <= p_cc ? :CC : :NC

            # Perform the interaction
            interact!(particle, interaction; rng=rng)

            # If CC interaction created a charged lepton, propagate it
            if interaction == :CC && particle.survived
                propagate_charged_lepton!(clp, particle, track)
            end

        else
            # Charged lepton: check if essentially at exit
            if 1.0 - particle.position < 1e-5
                particle.position = 1.0
                return particle
            end

            # Store decay position and decay
            particle.decay_position = particle.position
            decay!(particle; rng=rng)
        end
    end

    return particle
end

"""
    run_mc(energies, thetas, body, xs; kwargs...)

Run Monte Carlo simulation for multiple particles.

# Arguments
- `energies`: Vector of initial energies (eV)
- `thetas`: Vector of incident angles (radians)
- `body`: Body to propagate through
- `xs`: CrossSections object
- `flavor`: Initial neutrino flavor (default: NuTau)
- `secondaries`: Track secondary particles (default: true)
- `losses`: Include energy losses (default: true)
- `track_type`: :chord or :radial (default: :chord)
- `seed`: Random seed (default: nothing)

# Returns
Vector of result tuples: (E_initial, E_final, theta, n_cc, n_nc, final_id, position)
"""
function run_mc(
    energies::AbstractVector{<:Real},
    thetas::AbstractVector{<:Real},
    body::AbstractBody,
    xs::CrossSections;
    flavor::ParticleType=NuTau,
    secondaries::Bool=true,
    losses::Bool=true,
    track_type::Symbol=:chord,
    seed::Union{Nothing, Integer}=nothing,
    condition=nothing
)
    n_events = length(energies)
    if length(thetas) != n_events
        throw(ArgumentError("energies and thetas must have the same length"))
    end

    # Set up RNG (both Julia and PROPOSAL)
    rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)
    if !isnothing(seed)
        ChargedLeptonPropagationModule.seed_proposal!(seed)
    end

    # Create propagator
    clp = if body isa AbstractSphericalBody
        SphericalBodyPropagator(body)
    else
        SlabPropagator(body)
    end

    # Results storage
    results = Vector{NamedTuple{
        (:E_initial, :E_final, :theta, :n_cc, :n_nc, :id, :position),
        Tuple{Float64, Float64, Float64, Int, Int, Int, Float64}
    }}(undef, n_events)

    # Run simulation for each event
    for i in 1:n_events
        energy = Float64(energies[i])
        theta = Float64(thetas[i])

        # Create track
        track = if track_type == :chord
            Chord(theta)
        elseif track_type == :radial
            Radial()
        else
            error("Unknown track type: $track_type")
        end

        # Create and propagate particle
        particle = Particle(flavor, energy, 0.0, xs;
                           secondaries=secondaries, losses=losses)
        propagate!(particle, track, body, clp; condition=condition, rng=rng)

        # Store result
        results[i] = (
            E_initial = particle.initial_energy,
            E_final = particle.energy,
            theta = theta,
            n_cc = particle.n_cc,
            n_nc = particle.n_nc,
            id = Int(particle.id),
            position = particle.position
        )
    end

    return results
end

"""
    run_mc_parallel(energies, thetas, body, xs; nthreads=Threads.nthreads(), kwargs...)

Run Monte Carlo simulation in parallel using Julia threads.

Same arguments as `run_mc`, plus `nthreads` to control parallelism.
Each thread gets its own RNG seeded from the base seed.
"""
function run_mc_parallel(
    energies::AbstractVector{<:Real},
    thetas::AbstractVector{<:Real},
    body::AbstractBody,
    xs::CrossSections;
    nthreads::Int=Threads.nthreads(),
    flavor::ParticleType=NuTau,
    secondaries::Bool=true,
    losses::Bool=true,
    track_type::Symbol=:chord,
    seed::Union{Nothing, Integer}=nothing,
    condition=nothing
)
    n_events = length(energies)
    if length(thetas) != n_events
        throw(ArgumentError("energies and thetas must have the same length"))
    end

    # Results storage
    results = Vector{NamedTuple{
        (:E_initial, :E_final, :theta, :n_cc, :n_nc, :id, :position),
        Tuple{Float64, Float64, Float64, Int, Int, Int, Float64}
    }}(undef, n_events)

    # Base seed for reproducibility
    base_seed = isnothing(seed) ? rand(UInt64) : UInt64(seed)
    if !isnothing(seed)
        ChargedLeptonPropagationModule.seed_proposal!(seed)
    end

    Threads.@threads for i in 1:n_events
        # Thread-local RNG for reproducibility
        rng = MersenneTwister(base_seed + UInt64(i))

        energy = Float64(energies[i])
        theta = Float64(thetas[i])

        # Create track
        track = if track_type == :chord
            Chord(theta)
        elseif track_type == :radial
            Radial()
        else
            error("Unknown track type: $track_type")
        end

        # Create propagator (thread-local)
        clp = if body isa AbstractSphericalBody
            SphericalBodyPropagator(body)
        else
            SlabPropagator(body)
        end

        # Create and propagate particle
        particle = Particle(flavor, energy, 0.0, xs;
                           secondaries=secondaries, losses=losses)
        propagate!(particle, track, body, clp; condition=condition, rng=rng)

        # Store result
        results[i] = (
            E_initial = particle.initial_energy,
            E_final = particle.energy,
            theta = theta,
            n_cc = particle.n_cc,
            n_nc = particle.n_nc,
            id = Int(particle.id),
            position = particle.position
        )
    end

    return results
end

# Export main functions
export propagate!, run_mc, run_mc_parallel
