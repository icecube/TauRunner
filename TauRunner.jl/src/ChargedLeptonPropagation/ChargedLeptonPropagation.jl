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
using PROPOSAL

using ..TauRunner: AbstractChargedLeptonPropagator
using ..TauRunner: AbstractBody, AbstractSphericalBody, AbstractSlabBody
using ..TauRunner: AbstractTrack, AbstractSphericalTrack, AbstractSlabTrack
using ..TauRunner: ParticleType, is_charged_lepton
using ..TauRunner: Tau, AntiTau, Muon, AntiMuon, Electron, Positron
using ..TauRunner: NuTau, NuTauBar, NuMu, NuMuBar
using ..TauRunner: units, PhysicsConstants
using ..TauRunner.Bodies: get_density, radius, layer_boundaries
using ..TauRunner.Tracks: x_to_d, d_to_x, x_to_r

# Type-stable cache for PROPOSAL function references.
mutable struct ProposalFunctions
    # Hot-path functions (called per propagation step)
    ParticleState::Any
    propagate::Any
    get_final_state::Any
    get_energy::Any
    get_propagated_distance::Any
    # Propagator creation (called once per particle type, then cached)
    create_propagator_tauminus::Any
    create_propagator_tauplus::Any
    create_propagator_muminus::Any
    create_propagator_muplus::Any
    create_propagator_eminus::Any
    create_propagator_eplus::Any
    # Particle type constants
    PARTICLE_TYPE_TAUMINUS::Any
    PARTICLE_TYPE_TAUPLUS::Any
    PARTICLE_TYPE_MUMINUS::Any
    PARTICLE_TYPE_MUPLUS::Any
    PARTICLE_TYPE_EMINUS::Any
    PARTICLE_TYPE_EPLUS::Any
    # RNG seeding
    set_random_seed::Any
end

ProposalFunctions() = ProposalFunctions(ntuple(_ -> nothing, fieldcount(ProposalFunctions))...)

const PROPOSAL_FN_STRUCT = ProposalFunctions()

# Legacy Dict interface for backward compatibility (used in SphericalBodyPropagator.jl)
const PROPOSAL_FN = Dict{Symbol, Any}()

function _sync_proposal_fn_dict!()
    for name in fieldnames(ProposalFunctions)
        PROPOSAL_FN[name] = getfield(PROPOSAL_FN_STRUCT, name)
    end
end

function _populate_proposal_functions!()
    M = PROPOSAL
    PROPOSAL_FN_STRUCT.ParticleState = getfield(M, :ParticleState)
    PROPOSAL_FN_STRUCT.propagate = getfield(M, :propagate)
    PROPOSAL_FN_STRUCT.get_final_state = getfield(M, :get_final_state)
    PROPOSAL_FN_STRUCT.get_energy = getfield(M, :get_energy)
    PROPOSAL_FN_STRUCT.get_propagated_distance = getfield(M, :get_propagated_distance)

    PROPOSAL_FN_STRUCT.create_propagator_tauminus = getfield(M, :create_propagator_tauminus)
    PROPOSAL_FN_STRUCT.create_propagator_tauplus = getfield(M, :create_propagator_tauplus)
    PROPOSAL_FN_STRUCT.create_propagator_muminus = getfield(M, :create_propagator_muminus)
    PROPOSAL_FN_STRUCT.create_propagator_muplus = getfield(M, :create_propagator_muplus)
    PROPOSAL_FN_STRUCT.create_propagator_eminus = getfield(M, :create_propagator_eminus)
    PROPOSAL_FN_STRUCT.create_propagator_eplus = getfield(M, :create_propagator_eplus)

    PROPOSAL_FN_STRUCT.PARTICLE_TYPE_TAUMINUS = getfield(M, :PARTICLE_TYPE_TAUMINUS)
    PROPOSAL_FN_STRUCT.PARTICLE_TYPE_TAUPLUS = getfield(M, :PARTICLE_TYPE_TAUPLUS)
    PROPOSAL_FN_STRUCT.PARTICLE_TYPE_MUMINUS = getfield(M, :PARTICLE_TYPE_MUMINUS)
    PROPOSAL_FN_STRUCT.PARTICLE_TYPE_MUPLUS = getfield(M, :PARTICLE_TYPE_MUPLUS)
    PROPOSAL_FN_STRUCT.PARTICLE_TYPE_EMINUS = getfield(M, :PARTICLE_TYPE_EMINUS)
    PROPOSAL_FN_STRUCT.PARTICLE_TYPE_EPLUS = getfield(M, :PARTICLE_TYPE_EPLUS)
    PROPOSAL_FN_STRUCT.set_random_seed = getfield(M, :set_random_seed)

    _sync_proposal_fn_dict!()
end

function __init__()
    if !PROPOSAL.is_library_available()
        error("""
        PROPOSAL.jl native library not available.

        TauRunner.jl requires PROPOSAL for charged lepton propagation.
        For instructions on building and installing the PROPOSAL native library,
        see the PROPOSAL.jl README.
        """)
    end
    _populate_proposal_functions!()
    @info "PROPOSAL.jl loaded successfully"
end

export ChargedLeptonPropagator, SphericalBodyPropagator, SlabPropagator
export propagate_charged_lepton!, seed_proposal!

"""
    seed_proposal!(seed::Integer)

Seed PROPOSAL's internal random number generator for reproducibility.
"""
function seed_proposal!(seed::Integer)
    PROPOSAL.set_random_seed(Int(seed))
    return nothing
end

"""
    suppress_proposal_warnings(f)

Execute `f()` while suppressing PROPOSAL C++ warnings (written to stdout).
PROPOSAL's JSON config API produces occasional "Maximum number of iterations
exceeded in Bisection" warnings that are expected and harmless.
"""
function suppress_proposal_warnings(f)
    return f()
end

include("ConfigGeneration.jl")
include("AbstractPropagator.jl")
include("SphericalBodyPropagator.jl")
include("SlabPropagator.jl")

end # module ChargedLeptonPropagationModule
