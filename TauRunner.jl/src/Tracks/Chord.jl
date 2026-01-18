"""
Chord track through a spherical body.

A chord is a straight line path through a sphere, parameterized by the
incident angle θ (theta), where:
- θ = 0 corresponds to passing through the center (nadir)
- θ = π/2 corresponds to tangent to the surface (Earth-skimming)

The track parameter x ∈ [0, 1] represents progress along the chord,
where x=0 is the entry point and x=1 is the exit point.
"""

"""
    Chord{T} <: AbstractSphericalTrack

A chord trajectory through a spherical body.

# Fields
- `theta::T`: Incident angle in radians (0 = through center, π/2 = tangent)
- `depth::T`: Detector depth as fraction of radius (0 = surface, 1 = center)
- `cos_theta::T`: Precomputed cos(θ)
- `sin_theta::T`: Precomputed sin(θ)
- `l1::T`: Geometric parameter (entry to closest approach)
- `l2::T`: Geometric parameter (closest approach to exit)
- `total_length::T`: Total chord length in units of radius
"""
struct Chord{T<:Real} <: AbstractSphericalTrack
    theta::T
    depth::T
    cos_theta::T
    sin_theta::T
    l1::T
    l2::T
    total_length::T
end

"""
    Chord(theta; depth=0.0)

Construct a Chord track.

# Arguments
- `theta`: Incident angle in radians
- `depth`: Detector depth as fraction of radius (default: 0.0 = surface)

# Examples
```julia
# Through the center (nadir)
chord = Chord(0.0)

# 45 degree angle
chord = Chord(π/4)

# Earth-skimming at surface
chord = Chord(π/2 - 0.01)

# With detector at 1 km depth in a 6371 km Earth
chord = Chord(π/4, depth=1.0/6371.0)
```
"""
function Chord(theta::Real; depth::Real=0.0)
    T = Float64

    if !(0 <= theta <= π)
        @warn "Theta should be between 0 and π radians. Got $theta."
    end

    c = cos(theta)
    s = sin(theta)

    # Geometric calculations (see Python implementation)
    # l1 = distance from entry to point of closest approach (in units of radius)
    # l2 = distance from closest approach to exit
    l1 = sqrt(c^2 + 2*depth*s^2 - depth^2*s^2)
    l2 = (1 - depth) * c
    total = l1 + l2

    return Chord{T}(T(theta), T(depth), T(c), T(s), T(l1), T(l2), T(total))
end

# Coordinate transformation functions

"""
    x_to_r(chord::Chord, x::Real)

Convert track parameter x to normalized radius r.
"""
function x_to_r(chord::Chord{T}, x::Real) where T
    t = chord.total_length
    l1 = chord.l1
    return T(sqrt(1 - 2*x*l1*t + x^2*t^2))
end

"""
    x_to_r_prime(chord::Chord, x::Real)

Derivative of r with respect to x.
"""
function x_to_r_prime(chord::Chord{T}, x::Real) where T
    r = x_to_r(chord, x)
    t = chord.total_length
    l1 = chord.l1
    num = -l1*t + x*t^2

    if num == 0 && r == 0
        return zero(T)
    else
        return T(num / r)
    end
end

"""
    r_to_x(chord::Chord, r::Real)

Convert normalized radius r to track parameter x.
May return a tuple of two values if the track crosses the radius twice.
"""
function r_to_x(chord::Chord{T}, r::Real) where T
    t = chord.total_length
    l1 = chord.l1

    discriminant = r^2 - 1 + l1^2
    if discriminant < 0
        throw(ArgumentError("Radius $r is not reachable on this chord"))
    end

    sqrt_disc = sqrt(discriminant)
    val1 = (l1 - sqrt_disc) / t
    val2 = (l1 + sqrt_disc) / t

    # Return single value if only one intersection, or both if two
    if val2 > 1 || abs(val1 - val2) < 1e-10
        return T(val1)
    else
        return (T(val1), T(val2))
    end
end

"""
    x_to_d(chord::Chord, x::Real)

Convert track parameter x to distance traveled (in units of body radius).
"""
function x_to_d(chord::Chord{T}, x::Real) where T
    return T(chord.total_length * x)
end

"""
    d_to_x(chord::Chord, d::Real)

Convert distance traveled to track parameter x.
"""
function d_to_x(chord::Chord{T}, d::Real) where T
    return T(d / chord.total_length)
end

"""
    x_to_d_prime(chord::Chord, x::Real)

Derivative of distance with respect to x (constant for chord).
"""
function x_to_d_prime(chord::Chord{T}, x::Real) where T
    return chord.total_length
end

"""
    x_to_cartesian_direction(chord::Chord, x::Real)

Get the direction vector in Cartesian coordinates at position x.
Returns a 3-tuple (dx, dy, dz) representing the unit direction vector.
"""
function x_to_cartesian_direction(chord::Chord{T}, x::Real) where T
    # Direction is constant along a chord
    dx = -chord.sin_theta
    dy = zero(T)
    dz = chord.cos_theta
    norm = sqrt(dx^2 + dy^2 + dz^2)
    return (T(dx/norm), T(dy/norm), T(dz/norm))
end

"""
Integrand for column depth calculation along a chord.
"""
function integrand_column_depth(chord::Chord, body::AbstractSphericalBody)
    function f(x)
        r = x_to_r(chord, x)
        density = get_density(body, r)
        dr_dx = abs(x_to_r_prime(chord, x))
        return density * dr_dx * radius(body)
    end
    return f
end

function Base.hash(chord::Chord, h::UInt)
    return hash((chord.theta, chord.depth), h)
end

function Base.:(==)(a::Chord, b::Chord)
    return a.theta == b.theta && a.depth == b.depth
end

function Base.show(io::IO, chord::Chord)
    theta_deg = rad2deg(chord.theta)
    print(io, "Chord(θ=$(round(theta_deg, digits=2))°, depth=$(round(chord.depth, sigdigits=3)))")
end
