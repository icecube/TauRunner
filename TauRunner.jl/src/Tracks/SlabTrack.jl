"""
Slab track for flat geometry.

A straight-line track through a slab body.
The track parameter x ∈ [0, 1] represents progress through the slab.
"""

"""
    SlabTrack{T} <: AbstractSlabTrack

A straight-line trajectory through a slab body.

# Fields
- `theta::T`: Angle with respect to normal (0 = perpendicular)
- `depth::T`: Starting depth as fraction of total length
"""
struct SlabTrack{T<:Real} <: AbstractSlabTrack
    theta::T
    depth::T
    cos_theta::T
    path_length_factor::T  # 1/cos(theta), accounts for angled path
end

"""
    SlabTrack(; theta=0.0, depth=0.0)

Construct a SlabTrack.

# Arguments
- `theta`: Angle with respect to slab normal in radians (default: 0.0 = perpendicular)
- `depth`: Starting depth as fraction of total length (default: 0.0)
"""
function SlabTrack(; theta::Real=0.0, depth::Real=0.0)
    T = Float64
    c = cos(theta)
    path_factor = 1.0 / c  # Path length is longer for angled tracks
    return SlabTrack{T}(T(theta), T(depth), T(c), T(path_factor))
end

# Coordinate transformation functions
# For slab, x directly maps to position along the slab

"""
    x_to_d(track::SlabTrack, x::Real)

Convert track parameter x to distance traveled (in units of body length).
"""
function x_to_d(track::SlabTrack{T}, x::Real) where T
    return T(x * track.path_length_factor)
end

"""
    d_to_x(track::SlabTrack, d::Real)

Convert distance traveled to track parameter x.
"""
function d_to_x(track::SlabTrack{T}, d::Real) where T
    return T(d / track.path_length_factor)
end

"""
    x_to_d_prime(track::SlabTrack, x::Real)

Derivative of distance with respect to x.
"""
function x_to_d_prime(track::SlabTrack{T}, x::Real) where T
    return track.path_length_factor
end

"""
    x_to_cartesian_direction(track::SlabTrack, x::Real)

Get the direction vector in Cartesian coordinates.
"""
function x_to_cartesian_direction(track::SlabTrack{T}, x::Real) where T
    s = sin(track.theta)
    c = track.cos_theta
    return (T(s), zero(T), T(c))
end

"""
Integrand for column depth calculation along a slab track.
"""
function integrand_column_depth(track::SlabTrack, body::AbstractSlabBody)
    function f(x)
        density = get_density(body, x)
        return density * track.path_length_factor * length(body)
    end
    return f
end

function Base.hash(track::SlabTrack, h::UInt)
    return hash((track.theta, track.depth), h)
end

function Base.:(==)(a::SlabTrack, b::SlabTrack)
    return a.theta == b.theta && a.depth == b.depth
end

function Base.show(io::IO, track::SlabTrack)
    theta_deg = rad2deg(track.theta)
    print(io, "SlabTrack(θ=$(round(theta_deg, digits=2))°, depth=$(round(track.depth, sigdigits=3)))")
end
