"""
Radial track through a spherical body.

A radial track goes from the surface toward the center of the sphere.
The track parameter x ∈ [0, 1] represents progress inward, where:
- x = 0 is the surface (r = 1)
- x = 1 is the center (r = 0)
"""

"""
    Radial{T} <: AbstractSphericalTrack

A radial trajectory through a spherical body (surface to center).

# Fields
- `depth::T`: Starting depth as fraction of radius (usually 0)
"""
struct Radial{T<:Real} <: AbstractSphericalTrack
    depth::T
end

"""
    Radial(; depth=0.0)

Construct a Radial track.

# Arguments
- `depth`: Starting depth as fraction of radius (default: 0.0 = surface)
"""
function Radial(; depth::Real=0.0)
    return Radial{Float64}(Float64(depth))
end

# Coordinate transformation functions

"""
    x_to_r(track::Radial, x::Real)

Convert track parameter x to normalized radius r.
For radial track: r = 1 - x (moving from surface toward center).
"""
function x_to_r(track::Radial{T}, x::Real) where T
    return T(1 - track.depth - x * (1 - track.depth))
end

"""
    x_to_r_prime(track::Radial, x::Real)

Derivative of r with respect to x (constant for radial track).
"""
function x_to_r_prime(track::Radial{T}, x::Real) where T
    return T(-(1 - track.depth))
end

"""
    r_to_x(track::Radial, r::Real)

Convert normalized radius r to track parameter x.
"""
function r_to_x(track::Radial{T}, r::Real) where T
    if !(0 <= r <= 1 - track.depth)
        throw(ArgumentError("Radius $r is out of range for this radial track"))
    end
    return T((1 - track.depth - r) / (1 - track.depth))
end

"""
    x_to_d(track::Radial, x::Real)

Convert track parameter x to distance traveled (in units of body radius).
"""
function x_to_d(track::Radial{T}, x::Real) where T
    return T(x * (1 - track.depth))
end

"""
    d_to_x(track::Radial, d::Real)

Convert distance traveled to track parameter x.
"""
function d_to_x(track::Radial{T}, d::Real) where T
    return T(d / (1 - track.depth))
end

"""
    x_to_d_prime(track::Radial, x::Real)

Derivative of distance with respect to x (constant for radial track).
"""
function x_to_d_prime(track::Radial{T}, x::Real) where T
    return T(1 - track.depth)
end

"""
    x_to_cartesian_direction(track::Radial, x::Real)

Get the direction vector in Cartesian coordinates at position x.
Returns a 3-tuple (dx, dy, dz) pointing toward the center.
"""
function x_to_cartesian_direction(track::Radial{T}, x::Real) where T
    # Pointing straight down (toward center)
    return (zero(T), zero(T), T(-1))
end

"""
Integrand for column depth calculation along a radial track.
"""
function integrand_column_depth(track::Radial, body::AbstractSphericalBody)
    function f(x)
        r = x_to_r(track, x)
        density = get_density(body, r)
        # dr/dx is constant and negative, we want absolute value
        dr_dx = abs(x_to_r_prime(track, x))
        return density * dr_dx * radius(body)
    end
    return f
end

function Base.hash(track::Radial, h::UInt)
    return hash(track.depth, h)
end

function Base.:(==)(a::Radial, b::Radial)
    return a.depth == b.depth
end

function Base.show(io::IO, track::Radial)
    print(io, "Radial(depth=$(round(track.depth, sigdigits=3)))")
end
