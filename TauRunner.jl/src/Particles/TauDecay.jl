"""
Tau decay parameterization and secondary neutrino sampling.

The tau lepton decays with:
- ~17.8% to electron + neutrinos (leptonic)
- ~17.4% to muon + neutrinos (leptonic)
- ~64.8% to hadrons + tau neutrino (hadronic)

For leptonic decays, we track the secondary anti-neutrinos.
"""

# Tau decay branching ratios
const TAU_BR_ELECTRON = 0.18  # τ → e + νe + ντ
const TAU_BR_MUON = 0.18      # τ → μ + νμ + ντ (total leptonic ~0.36)
const TAU_BR_HADRONIC = 0.64  # τ → hadrons + ντ

# Energy fractions retained by the tau neutrino in hadronic decays
# These are sampled from the decay distribution
const TAU_DECAY_FRACTIONS = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
const TAU_DECAY_WEIGHTS = [0.15, 0.15, 0.12, 0.10, 0.10, 0.10, 0.08, 0.08, 0.07, 0.05]

# TODO: Load actual CDF splines for secondary neutrino energy fractions
# For now, use simple approximations

"""
    sample_secondary_energy_fraction(rng::AbstractRNG, decay_type::Symbol)

Sample the energy fraction for secondary neutrinos from leptonic tau decays.

# Arguments
- `rng`: Random number generator
- `decay_type`: :electron or :muon

# Returns
Energy fraction for the secondary (anti)neutrino
"""
function sample_secondary_energy_fraction(rng::AbstractRNG, decay_type::Symbol)
    # Simplified sampling - in full implementation, use CDF splines
    # The secondary neutrino typically carries a fraction of the tau energy
    # distributed according to the 3-body decay kinematics

    # Simple approximation: uniform in [0.1, 0.5]
    return 0.1 + 0.4 * rand(rng)
end

"""
    decay!(particle::Particle; rng=Random.default_rng())

Perform a tau or muon decay, updating particle state.

For tau decays:
- Sample the decay channel (leptonic vs hadronic)
- For leptonic decays, add secondary neutrino to basket
- Convert tau to tau neutrino with sampled energy fraction

For muon/electron:
- Mark particle as not survived (absorbed)
"""
function decay!(particle::Particle{T}; rng::AbstractRNG=Random.default_rng()) where T
    abs_id = abs(Int(particle.id))
    sign_id = sign(Int(particle.id))

    if abs_id in (12, 14, 16)
        error("Neutrinos don't decay in this simulation")
    end

    if abs_id == 15  # Tau
        particle.n_decay += 1

        if particle.include_secondaries
            # Sample decay channel
            p0 = rand(rng)

            if p0 < TAU_BR_ELECTRON
                # τ → e + ν̄e + ντ (or τ̄ → ē + νe + ν̄τ)
                # Sample energy of secondary anti-nue (for τ-) or nue (for τ+)
                frac = sample_secondary_energy_fraction(rng, :electron)
                E_secondary = frac * particle.energy
                secondary_id = ParticleType(-sign_id * 12)  # opposite sign from tau

                push!(particle.basket, SecondaryParticle{T}(
                    secondary_id,
                    T(E_secondary),
                    particle.position
                ))

            elseif p0 < TAU_BR_ELECTRON + TAU_BR_MUON
                # τ → μ + ν̄μ + ντ (or τ̄ → μ̄ + νμ + ν̄τ)
                frac = sample_secondary_energy_fraction(rng, :muon)
                E_secondary = frac * particle.energy
                secondary_id = ParticleType(-sign_id * 14)  # opposite sign from tau

                push!(particle.basket, SecondaryParticle{T}(
                    secondary_id,
                    T(E_secondary),
                    particle.position
                ))
            end
            # else: hadronic decay, no tracked secondary
        end

        # Sample energy fraction retained by tau neutrino
        z = sample_weighted(rng, TAU_DECAY_FRACTIONS, TAU_DECAY_WEIGHTS)
        particle.energy = T(z * particle.energy)

        # Convert to tau neutrino
        particle.id = ParticleType(sign_id * 16)
        particle.mass, particle.lifetime = get_particle_properties(particle.id)

    elseif abs_id in (11, 13)  # Electron or Muon
        # These effectively stop/are absorbed
        particle.survived = false
    else
        error("Unknown particle type for decay: $(particle.id)")
    end

    return nothing
end
