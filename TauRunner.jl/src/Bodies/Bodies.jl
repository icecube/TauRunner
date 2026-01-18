"""
Bodies module for TauRunner.jl

Defines physical bodies through which particles propagate (Earth, Sun, custom spheres, slabs).
"""
module Bodies

using QuadGK
using ..TauRunner: AbstractBody, AbstractSphericalBody, AbstractSlabBody
using ..TauRunner: units, PhysicsConstants

export SphericalBody, LayeredConstantSlab
export construct_earth, construct_sun, construct_uniform_earth
export get_density, get_average_density, radius, radius_km, layer_boundaries

include("SphericalBody.jl")
include("Earth.jl")
include("Slab.jl")

end # module Bodies
