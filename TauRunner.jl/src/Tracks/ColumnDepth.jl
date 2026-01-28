"""
Column depth calculation utilities for TauRunner.jl

Provides functions for computing and caching column depth along tracks.
"""

# LRUCache is already imported by the parent module

# Cache for column depth splines
# Key: hash of (track, body), Value: (total_depth, x_to_X_spline, X_to_x_spline)
const COLUMN_DEPTH_CACHE = LRU{UInt64, Tuple{Float64, Any, Any}}(maxsize=100)

"""
    compute_column_depth_hash(track::AbstractTrack, body::AbstractBody)

Compute a hash key for caching column depth calculations.
"""
function compute_column_depth_hash(track::AbstractTrack, body::AbstractBody)
    return hash((typeof(track), track, typeof(body), body))
end

"""
    integrand_column_depth(track::AbstractTrack, body::AbstractBody)

Return a function that computes the column depth integrand at position x.
The integrand is: ρ(r(x)) * |dr/dx| * R
where R is the body radius/length.
"""
function integrand_column_depth end

"""
    compute_column_depth_spline(track::AbstractTrack, body::AbstractBody; n_points=1000)

Compute the column depth spline for a given track through a body.

Returns a tuple of:
1. Total column depth
2. Interpolation function x_to_X (position to accumulated depth)
3. Interpolation function X_to_x (accumulated depth to position)
"""
function compute_column_depth_spline(
    track::AbstractTrack,
    body::AbstractBody;
    n_points::Int=1000
)
    # Sample points along the track
    xs = range(0.0, 1.0, length=n_points)
    Xs = zeros(Float64, n_points)

    # Get the integrand function
    f = integrand_column_depth(track, body)

    # Compute cumulative column depth at each point
    for i in 2:n_points
        integral, _ = quadgk(f, xs[i-1], xs[i])
        Xs[i] = Xs[i-1] + integral
    end

    total_depth = Xs[end]

    # Create interpolation for x -> X
    x_to_X_interp = linear_interpolation(xs, Xs)

    # Create interpolation for X -> x (need to handle monotonicity)
    # Filter out any non-strictly-increasing points
    valid_indices = [1]
    for i in 2:n_points
        if Xs[i] > Xs[valid_indices[end]]
            push!(valid_indices, i)
        end
    end

    if length(valid_indices) >= 2
        X_to_x_interp = linear_interpolation(Xs[valid_indices], xs[valid_indices])
    else
        error("Column depth is degenerate (not strictly increasing) for this track/body combination. " *
              "Total column depth: $total_depth. Check body density and track geometry.")
    end

    return (total_depth, x_to_X_interp, X_to_x_interp)
end

"""
    get_column_depth_data(track::AbstractTrack, body::AbstractBody)

Get cached column depth data, computing if necessary.
"""
function get_column_depth_data(track::AbstractTrack, body::AbstractBody)
    key = compute_column_depth_hash(track, body)

    return get!(COLUMN_DEPTH_CACHE, key) do
        compute_column_depth_spline(track, body)
    end
end

"""
    total_column_depth(track::AbstractTrack, body::AbstractBody)

Compute the total column depth along the entire track through the body.

# Returns
Total column depth in natural units (eV^3, since it's density × length).
"""
function total_column_depth(track::AbstractTrack, body::AbstractBody)
    data = get_column_depth_data(track, body)
    return data[1]
end

"""
    x_to_X(track::AbstractTrack, body::AbstractBody, x::Real)

Convert track parameter x to accumulated column depth X.

# Arguments
- `track`: Track object
- `body`: Body object
- `x`: Track parameter (0 = entry, 1 = exit)

# Returns
Accumulated column depth from entry to position x.
"""
function x_to_X(track::AbstractTrack, body::AbstractBody, x::Real)
    if !(0 <= x <= 1)
        throw(ArgumentError("x must be between 0 and 1, got $x"))
    end
    data = get_column_depth_data(track, body)
    return data[2](x)
end

"""
    X_to_x(track::AbstractTrack, body::AbstractBody, X::Real)

Convert accumulated column depth X to track parameter x.

# Arguments
- `track`: Track object
- `body`: Body object
- `X`: Accumulated column depth

# Returns
Track parameter x at which the accumulated depth equals X.

# Throws
`ArgumentError` if X exceeds the total column depth.
"""
function X_to_x(track::AbstractTrack, body::AbstractBody, X::Real)
    data = get_column_depth_data(track, body)
    total = data[1]

    if X > total
        throw(ArgumentError("Column depth $X exceeds total $total"))
    end

    return data[3](X)
end

"""
    clear_column_depth_cache!()

Clear the column depth cache. Useful for testing or memory management.
"""
function clear_column_depth_cache!()
    empty!(COLUMN_DEPTH_CACHE)
    return nothing
end
