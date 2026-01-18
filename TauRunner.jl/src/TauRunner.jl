"""
TauRunner.jl - Monte Carlo simulation of ultra-high-energy neutrino propagation

A Julia port of the TauRunner Python package for propagating neutrinos through
matter (Earth, Sun, arbitrary spherical bodies) including interactions and energy losses.

Authors: Ibrahim Safa, Carlos A. Arguelles, Jeff Lazar, Alex Pizzuto
"""
module TauRunner

using LinearAlgebra
using Random
using StaticArrays

# Include constants and units first (no dependencies)
include("Constants.jl")
include("Units.jl")

# Abstract type hierarchy
# Bodies
abstract type AbstractBody end
abstract type AbstractSphericalBody <: AbstractBody end
abstract type AbstractSlabBody <: AbstractBody end

# Tracks
abstract type AbstractTrack end
abstract type AbstractSphericalTrack <: AbstractTrack end
abstract type AbstractSlabTrack <: AbstractTrack end

# Cross-sections
abstract type AbstractCrossSections end

# Charged lepton propagators
abstract type AbstractChargedLeptonPropagator end

# Particle types (PDG IDs)
@enum ParticleType begin
    NuE = 12
    NuMu = 14
    NuTau = 16
    NuEBar = -12
    NuMuBar = -14
    NuTauBar = -16
    Electron = 11
    Muon = 13
    Tau = 15
    Positron = -11
    AntiMuon = -13
    AntiTau = -15
end

# Helper functions for particle types
const NEUTRINO_IDS = (NuE, NuMu, NuTau, NuEBar, NuMuBar, NuTauBar)
const CHARGED_LEPTON_IDS = (Electron, Muon, Tau, Positron, AntiMuon, AntiTau)
const LEPTON_IDS = (NEUTRINO_IDS..., CHARGED_LEPTON_IDS...)

is_neutrino(p::ParticleType) = p in NEUTRINO_IDS
is_charged_lepton(p::ParticleType) = p in CHARGED_LEPTON_IDS
is_antineutrino(p::ParticleType) = Int(p) < 0 && is_neutrino(p)

function charged_partner(p::ParticleType)
    id = Int(p)
    if abs(id) in (12, 14, 16)
        return ParticleType(sign(id) * (abs(id) - 1))
    else
        error("No charged partner for particle type $p")
    end
end

function neutrino_partner(p::ParticleType)
    id = Int(p)
    if abs(id) in (11, 13, 15)
        return ParticleType(sign(id) * (abs(id) + 1))
    else
        error("No neutrino partner for particle type $p")
    end
end

# Include submodules
include("Bodies/Bodies.jl")
include("Tracks/Tracks.jl")
include("CrossSections/CrossSections.jl")
include("Particles/Particles.jl")
include("ChargedLeptonPropagation/ChargedLeptonPropagation.jl")
include("Casino.jl")

# Re-export from submodules
using .Bodies
using .Tracks
using .CrossSectionsModule
using .Particles
using .ChargedLeptonPropagationModule

# Export main types and functions
export AbstractBody, AbstractSphericalBody, AbstractSlabBody
export AbstractTrack, AbstractSphericalTrack, AbstractSlabTrack
export AbstractCrossSections
export AbstractChargedLeptonPropagator
export ParticleType, NEUTRINO_IDS, CHARGED_LEPTON_IDS, LEPTON_IDS
export NuE, NuMu, NuTau, NuEBar, NuMuBar, NuTauBar
export Electron, Muon, Tau, Positron, AntiMuon, AntiTau
export is_neutrino, is_charged_lepton, is_antineutrino
export charged_partner, neutrino_partner

# Export units and constants
export units, PhysicsConstants

# Export from Bodies
export SphericalBody, construct_earth, get_density, get_average_density

# Export from Tracks
export Chord, Radial, SlabTrack
export x_to_r, r_to_x, x_to_d, d_to_x, total_column_depth, X_to_x, x_to_X

# Export from CrossSections
export CrossSections, XSModel, DIPOLE, CSMS
export total_cross_section, differential_cross_section

# Export from Particles
export Particle, SecondaryParticle
export get_proposed_depth_step, get_interaction_depth, interact!, decay!

# Export from ChargedLeptonPropagation
export ChargedLeptonPropagator, propagate_charged_lepton!

# Export main propagation function
export propagate!

end # module TauRunner
