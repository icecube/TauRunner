"""
CrossSections module for TauRunner.jl

Handles neutrino cross-section calculations for CC and NC interactions.
"""
module CrossSectionsModule

using Interpolations
using JLD2
using JSON3
using NPZ
using ..TauRunner: AbstractCrossSections
using ..TauRunner: units, PhysicsConstants

export CrossSections, XSModel, DIPOLE, CSMS
export total_cross_section, differential_cross_section

include("XSModel.jl")
include("Interpolation.jl")

end # module CrossSectionsModule
