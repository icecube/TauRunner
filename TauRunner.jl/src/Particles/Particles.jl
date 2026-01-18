"""
Particles module for TauRunner.jl

Defines particle state and interaction methods.
"""
module Particles

using Random
using ..TauRunner: ParticleType, AbstractCrossSections, AbstractBody, AbstractTrack
using ..TauRunner: is_neutrino, is_charged_lepton, charged_partner, neutrino_partner
using ..TauRunner: NuE, NuMu, NuTau, NuEBar, NuMuBar, NuTauBar
using ..TauRunner: Electron, Muon, Tau, Positron, AntiMuon, AntiTau
using ..TauRunner: units, PhysicsConstants
using ..TauRunner.CrossSectionsModule: total_cross_section, differential_cross_section

export Particle, SecondaryParticle
export get_proposed_depth_step, get_interaction_depth, get_total_interaction_depth
export interact!, decay!, reset!

include("Particle.jl")
include("TauDecay.jl")

end # module Particles
