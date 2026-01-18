"""
ChargedLeptonPropagation module for TauRunner.jl

Handles propagation of charged leptons (tau, muon) through matter,
including energy losses from ionization, bremsstrahlung, pair production,
and photonuclear interactions.

This module interfaces with PROPOSAL.jl for charged lepton propagation.
"""
module ChargedLeptonPropagationModule

using JSON3
using Random

using ..TauRunner: AbstractChargedLeptonPropagator
using ..TauRunner: AbstractBody, AbstractSphericalBody, AbstractSlabBody
using ..TauRunner: AbstractTrack, AbstractSphericalTrack, AbstractSlabTrack
using ..TauRunner: ParticleType, is_charged_lepton
using ..TauRunner: Tau, AntiTau, Muon, AntiMuon, Electron, Positron
using ..TauRunner: NuTau, NuTauBar, NuMu, NuMuBar
using ..TauRunner: units, PhysicsConstants
using ..TauRunner.Bodies: get_density, radius, layer_boundaries
using ..TauRunner.Tracks: x_to_d, d_to_x, x_to_r

# Try to import PROPOSAL.jl - it may not be available
const PROPOSAL_AVAILABLE = Ref(false)
const SIMPLIFIED_WARNING_SHOWN = Ref(false)

function __init__()
    try
        @eval using PROPOSAL
        PROPOSAL_AVAILABLE[] = true
        @info "PROPOSAL.jl loaded successfully"
    catch e
        PROPOSAL_AVAILABLE[] = false
        @warn """
        PROPOSAL.jl not available. Falling back to simplified charged lepton propagation.

        ⚠️  WARNING: The simplified propagation model is NOT suitable for physics analysis!

        Limitations of simplified model:
        • Uses average continuous energy loss (~2 MeV/(g/cm²)) instead of stochastic sampling
        • No bremsstrahlung, pair production, or photonuclear interaction sampling
        • No multiple scattering or angular deflection
        • Decay uses mean decay length instead of proper Monte Carlo sampling

        For accurate physics results, install PROPOSAL.jl:
            using Pkg; Pkg.add("PROPOSAL")
        """ exception=e
    end
end

"""
    warn_simplified_propagation()

Issue a one-time warning that simplified propagation is being used.
"""
function warn_simplified_propagation()
    if !SIMPLIFIED_WARNING_SHOWN[]
        SIMPLIFIED_WARNING_SHOWN[] = true
        @warn """
        Using simplified charged lepton propagation (PROPOSAL.jl not available).
        Results are NOT suitable for physics analysis. See TauRunner documentation for details.
        """
    end
end

export ChargedLeptonPropagator, SphericalBodyPropagator, SlabPropagator
export propagate_charged_lepton!, is_proposal_available

"""Check if PROPOSAL.jl is available."""
is_proposal_available() = PROPOSAL_AVAILABLE[]

include("ConfigGeneration.jl")
include("AbstractPropagator.jl")
include("SphericalBodyPropagator.jl")
include("SlabPropagator.jl")

end # module ChargedLeptonPropagationModule
