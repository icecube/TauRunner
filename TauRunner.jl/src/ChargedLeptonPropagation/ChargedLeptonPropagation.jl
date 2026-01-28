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

# Type-stable cache for PROPOSAL function references.
# Using a mutable struct with concrete `Any` fields avoids Dict lookup overhead
# on every call in hot loops, while keeping the late-binding semantics needed
# because PROPOSAL types are only known at runtime.
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
# Maps symbol keys to field access on the struct.
const PROPOSAL_FN = Dict{Symbol, Any}()

function _sync_proposal_fn_dict!()
    for name in fieldnames(ProposalFunctions)
        PROPOSAL_FN[name] = getfield(PROPOSAL_FN_STRUCT, name)
    end
end

function __init__()
    try
        @eval using PROPOSAL
        PROPOSAL_AVAILABLE[] = true

        # Cache function references in type-stable struct
        PROPOSAL_FN_STRUCT.ParticleState = @eval PROPOSAL.ParticleState
        PROPOSAL_FN_STRUCT.propagate = @eval PROPOSAL.propagate
        PROPOSAL_FN_STRUCT.get_final_state = @eval PROPOSAL.get_final_state
        PROPOSAL_FN_STRUCT.get_energy = @eval PROPOSAL.get_energy
        PROPOSAL_FN_STRUCT.get_propagated_distance = @eval PROPOSAL.get_propagated_distance

        PROPOSAL_FN_STRUCT.create_propagator_tauminus = @eval PROPOSAL.create_propagator_tauminus
        PROPOSAL_FN_STRUCT.create_propagator_tauplus = @eval PROPOSAL.create_propagator_tauplus
        PROPOSAL_FN_STRUCT.create_propagator_muminus = @eval PROPOSAL.create_propagator_muminus
        PROPOSAL_FN_STRUCT.create_propagator_muplus = @eval PROPOSAL.create_propagator_muplus
        PROPOSAL_FN_STRUCT.create_propagator_eminus = @eval PROPOSAL.create_propagator_eminus
        PROPOSAL_FN_STRUCT.create_propagator_eplus = @eval PROPOSAL.create_propagator_eplus

        PROPOSAL_FN_STRUCT.PARTICLE_TYPE_TAUMINUS = @eval PROPOSAL.PARTICLE_TYPE_TAUMINUS
        PROPOSAL_FN_STRUCT.PARTICLE_TYPE_TAUPLUS = @eval PROPOSAL.PARTICLE_TYPE_TAUPLUS
        PROPOSAL_FN_STRUCT.PARTICLE_TYPE_MUMINUS = @eval PROPOSAL.PARTICLE_TYPE_MUMINUS
        PROPOSAL_FN_STRUCT.PARTICLE_TYPE_MUPLUS = @eval PROPOSAL.PARTICLE_TYPE_MUPLUS
        PROPOSAL_FN_STRUCT.PARTICLE_TYPE_EMINUS = @eval PROPOSAL.PARTICLE_TYPE_EMINUS
        PROPOSAL_FN_STRUCT.PARTICLE_TYPE_EPLUS = @eval PROPOSAL.PARTICLE_TYPE_EPLUS
        PROPOSAL_FN_STRUCT.set_random_seed = @eval PROPOSAL.set_random_seed

        # Sync to Dict for code that still uses PROPOSAL_FN[key] syntax
        _sync_proposal_fn_dict!()

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
export propagate_charged_lepton!, is_proposal_available, seed_proposal!

"""Check if PROPOSAL.jl is available."""
is_proposal_available() = PROPOSAL_AVAILABLE[]

"""
    seed_proposal!(seed::Integer)

Seed PROPOSAL's internal random number generator for reproducibility.
"""
function seed_proposal!(seed::Integer)
    if PROPOSAL_AVAILABLE[]
        Base.invokelatest(PROPOSAL_FN_STRUCT.set_random_seed, Int(seed))
    end
    return nothing
end

"""
    suppress_proposal_warnings(f)

Execute `f()` while suppressing PROPOSAL C++ warnings (written to stdout).
PROPOSAL's JSON config API produces occasional "Maximum number of iterations
exceeded in Bisection" warnings that are expected and harmless.
"""
function suppress_proposal_warnings(f)
    old_stdout = ccall(:dup, Cint, (Cint,), 1)
    old_stderr = ccall(:dup, Cint, (Cint,), 2)
    devnull_fd = ccall(:open, Cint, (Cstring, Cint), "/dev/null", 1)
    ccall(:dup2, Cint, (Cint, Cint), devnull_fd, 1)
    ccall(:dup2, Cint, (Cint, Cint), devnull_fd, 2)
    ccall(:close, Cint, (Cint,), devnull_fd)
    try
        return f()
    finally
        ccall(:dup2, Cint, (Cint, Cint), old_stdout, 1)
        ccall(:dup2, Cint, (Cint, Cint), old_stderr, 2)
        ccall(:close, Cint, (Cint,), old_stdout)
        ccall(:close, Cint, (Cint,), old_stderr)
    end
end

include("ConfigGeneration.jl")
include("AbstractPropagator.jl")
include("SphericalBodyPropagator.jl")
include("SlabPropagator.jl")

end # module ChargedLeptonPropagationModule
