"""
Tau decay parameterization and secondary neutrino sampling.

Physics-based tau decay spectrum matching Python TauRunner implementation.
Based on the tau decay kinematics for different channels:
- Leptonic: τ → l + νl + ντ (3-body decay)
- Pion: τ → π + ντ (2-body decay)
- Rho: τ → ρ + ντ (2-body decay)
- A1: τ → a1 + ντ (2-body decay)
- Other hadronic: simplified flat distribution

References:
- Particle Data Group tau branching ratios
- Standard tau polarization effects (P = -1 for left-handed taus)
"""

# =============================================================================
# Mass ratios (m_meson / m_tau)^2
# =============================================================================
const R_PION = 0.07856^2   # (m_π / m_τ)²
const R_RHO = 0.43335^2    # (m_ρ / m_τ)²
const R_A1 = 0.70913^2     # (m_a1 / m_τ)²

# =============================================================================
# Branching ratios
# =============================================================================
const BR_LEPTON = 0.18     # τ → e/μ + neutrinos (each channel)
const BR_PION = 0.12       # τ → π + ντ
const BR_RHO = 0.26        # τ → ρ + ντ
const BR_A1 = 0.13         # τ → a1 + ντ
const BR_HAD = 0.13        # τ → other hadrons + ντ

# Total leptonic branching ratio (electron + muon)
const TAU_BR_ELECTRON = 0.18
const TAU_BR_MUON = 0.18
const TAU_BR_HADRONIC = 1.0 - TAU_BR_ELECTRON - TAU_BR_MUON

# =============================================================================
# Decay spectrum functions
# Each returns dN/dz where z = E_ντ / E_τ
# P is the tau polarization (-1 for left-handed taus from CC interactions)
# =============================================================================

"""
Leptonic decay spectrum: τ → l + νl + ντ (3-body decay)
"""
function tau_decay_to_lepton(z::Real, P::Real=-1.0)
    g0 = (5.0/3.0) - 3.0*z^2 + (4.0/3.0)*z^3
    g1 = (1.0/3.0) - 3.0*z^2 + (8.0/3.0)*z^3
    return g0 + P*g1
end

"""
Pion decay spectrum: τ → π + ντ (2-body decay)
"""
function tau_decay_to_pion(z::Real, P::Real=-1.0)
    if (1.0 - R_PION - z) > 0.0
        g0 = 1.0 / (1.0 - R_PION)
        g1 = -(2.0*z - 1.0 + R_PION) / (1.0 - R_PION)^2
        return g0 + P*g1
    else
        return 0.0
    end
end

"""
Rho decay spectrum: τ → ρ + ντ (2-body decay)
"""
function tau_decay_to_rho(z::Real, P::Real=-1.0)
    if (1.0 - R_RHO - z) > 0.0
        g0 = 1.0 / (1.0 - R_RHO)
        g1 = -((2.0*z - 1.0 + R_RHO) / (1.0 - R_RHO)) * ((1.0 - 2.0*R_RHO) / (1.0 + 2.0*R_RHO))
        return g0 + P*g1
    else
        return 0.0
    end
end

"""
A1 decay spectrum: τ → a1 + ντ (2-body decay)
"""
function tau_decay_to_a1(z::Real, P::Real=-1.0)
    if (1.0 - R_A1 - z) > 0.0
        g0 = 1.0 / (1.0 - R_A1)
        g1 = -((2.0*z - 1.0 + R_A1) / (1.0 - R_A1)) * ((1.0 - 2.0*R_A1) / (1.0 + 2.0*R_A1))
        return g0 + P*g1
    else
        return 0.0
    end
end

"""
Other hadronic decay spectrum (simplified flat distribution for z < 0.3)
"""
function tau_decay_to_hadrons(z::Real, P::Real=-1.0)
    if (0.3 - z) > 0.0
        return 1.0 / 0.3
    else
        return 0.0
    end
end

"""
Combined tau decay spectrum dN/dz for tau neutrino energy fraction z.
Includes all decay channels weighted by branching ratios.
"""
function tau_decay_spectrum(z::Real, P::Real=-1.0)
    spectrum = 0.0
    spectrum += 2.0 * BR_LEPTON * tau_decay_to_lepton(z, P)  # factor 2 for e + μ
    spectrum += BR_PION * tau_decay_to_pion(z, P)
    spectrum += BR_RHO * tau_decay_to_rho(z, P)
    spectrum += BR_A1 * tau_decay_to_a1(z, P)
    spectrum += BR_HAD * tau_decay_to_hadrons(z, P)
    return spectrum
end

# =============================================================================
# Precomputed decay distribution for efficient sampling
# Matches Python: 998 bins from 0.001 to 0.999
# =============================================================================
const N_DECAY_BINS = 998
const TAU_DECAY_FRACTIONS = collect(range(0.001, 0.999, length=N_DECAY_BINS))

# Compute weights from decay spectrum (P = -1 for left-handed taus)
const TAU_DECAY_WEIGHTS = let
    weights = [tau_decay_spectrum(z, -1.0) for z in TAU_DECAY_FRACTIONS]
    weights ./ sum(weights)  # Normalize
end

# Precompute CDF for efficient sampling
const TAU_DECAY_CDF = cumsum(TAU_DECAY_WEIGHTS)

"""
    sample_tau_decay_fraction(rng::AbstractRNG)

Sample the energy fraction z retained by the tau neutrino in tau decay.
Uses the physics-based decay spectrum with proper kinematics.

# Returns
Energy fraction z ∈ (0, 1) for the outgoing tau neutrino.
"""
function sample_tau_decay_fraction(rng::AbstractRNG)
    u = rand(rng)
    idx = searchsortedfirst(TAU_DECAY_CDF, u)
    idx = clamp(idx, 1, N_DECAY_BINS)
    return TAU_DECAY_FRACTIONS[idx]
end

# =============================================================================
# Secondary neutrino sampling (for leptonic decays)
# =============================================================================

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
    # TODO: Load actual CDF splines for secondary neutrino energy fractions
    # For now, use simple approximation: uniform in [0.1, 0.5]
    return 0.1 + 0.4 * rand(rng)
end

# =============================================================================
# Main decay function
# =============================================================================

"""
    decay!(particle::Particle; rng=Random.default_rng())

Perform a tau or muon decay, updating particle state.

For tau decays:
- Sample the decay channel (leptonic vs hadronic)
- For leptonic decays, add secondary neutrino to basket
- Convert tau to tau neutrino with energy sampled from physics-based spectrum

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

        # Sample energy fraction retained by tau neutrino using physics-based spectrum
        z = sample_tau_decay_fraction(rng)
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
