"""
Abstract interface for charged lepton propagators.

This defines the interface that concrete propagators (using PROPOSAL.jl) must implement.
"""

"""
    ChargedLeptonPropagator

Abstract type for charged lepton propagators.

Concrete subtypes must implement:
- `propagate_charged_lepton!(propagator, particle, track)`
"""
abstract type ChargedLeptonPropagator <: AbstractChargedLeptonPropagator end

"""
    propagate_charged_lepton!(propagator, particle, track)

Propagate a charged lepton through matter, updating its energy and position.

# Arguments
- `propagator`: A ChargedLeptonPropagator instance
- `particle`: The particle to propagate (modified in place)
- `track`: The track geometry

This function should:
1. Propagate the particle from its current position toward the exit
2. Handle energy losses (ionization, bremsstrahlung, etc.) via PROPOSAL
3. Check for decay and update particle state accordingly
4. Update particle.position and particle.energy

The particle may decay during propagation, in which case the decay products
should be handled appropriately.
"""
function propagate_charged_lepton! end

"""
    create_propagator(body::AbstractBody)

Factory function to create the appropriate propagator for a body type.
"""
function create_propagator(body::AbstractSphericalBody)
    return SphericalBodyPropagator(body)
end

function create_propagator(body::AbstractSlabBody)
    return SlabPropagator(body)
end
