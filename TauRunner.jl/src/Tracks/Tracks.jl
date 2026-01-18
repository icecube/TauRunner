"""
Tracks module for TauRunner.jl

Defines particle trajectories through bodies (chords through spheres, radial paths, slab tracks).
"""
module Tracks

using QuadGK
using Roots
using Interpolations
using LRUCache
using ..TauRunner: AbstractTrack, AbstractSphericalTrack, AbstractSlabTrack
using ..TauRunner: AbstractBody, AbstractSphericalBody, AbstractSlabBody
using ..TauRunner: units
using ..TauRunner.Bodies: get_density, radius, layer_boundaries

export Chord, Radial, SlabTrack
export x_to_r, r_to_x, x_to_d, d_to_x, x_to_d_prime, x_to_r_prime
export total_column_depth, X_to_x, x_to_X
export x_to_cartesian_direction

include("ColumnDepth.jl")
include("Chord.jl")
include("Radial.jl")
include("SlabTrack.jl")

end # module Tracks
