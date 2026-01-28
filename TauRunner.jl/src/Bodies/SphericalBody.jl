"""
SphericalBody type for representing layered spherical objects.

The density is specified as a function of normalized radius r ∈ [0, 1],
where r=0 is the center and r=1 is the surface.
"""

"""
    DensityFunction{T}

Wrapper that makes both constants and functions callable as density(r).
Analogous to the Python Callable class.
"""
struct DensityFunction{T<:Real, F}
    func::F
    scale::T
    is_callable::Bool
end

function DensityFunction(obj::F, scale::T=one(T)) where {F<:Function, T<:Real}
    DensityFunction{T, F}(obj, scale, true)
end

function DensityFunction(obj::T, scale::S=one(S)) where {T<:Real, S<:Real}
    R = promote_type(T, S)
    DensityFunction{R, T}(obj, R(scale), false)
end

function (df::DensityFunction{T})(r::Real) where T
    if df.is_callable
        return T(df.func(r) * df.scale)
    else
        return T(df.func * df.scale)
    end
end

"""
    SphericalBody{T, N} <: AbstractSphericalBody

A spherical body with N layers, each with its own density function.

# Fields
- `radius_natural::T`: Radius in natural units (eV^-1)
- `layer_boundaries::NTuple{N, T}`: Normalized boundaries (0 to 1) for each layer
- `density_functions::Vector{DensityFunction{T}}`: Density function for each layer
- `average_densities::Vector{T}`: Precomputed average density in each layer
- `name::Symbol`: Name identifier for the body
"""
struct SphericalBody{T<:Real} <: AbstractSphericalBody
    radius_natural::T
    layer_boundaries::Vector{T}
    density_functions::Vector{DensityFunction{T, <:Any}}
    average_densities::Vector{T}
    name::Symbol
end

"""
    SphericalBody(density, radius_km; layer_boundaries=nothing, name=:unnamed)

Construct a SphericalBody from density specification.

# Arguments
- `density`: Either a single density value/function, or a vector of (density, boundary) tuples
- `radius_km`: Radius of the body in kilometers

# Examples
```julia
# Uniform density sphere
body = SphericalBody(2.5, 6371.0, name=:uniform_earth)

# Layered sphere with density functions
layers = [
    (r -> 13.0 - 8.8*r, 0.19),  # Inner core
    (r -> 12.5 - 3.6*r, 0.55),  # Outer core
    (5.5, 1.0)                   # Mantle (constant)
]
body = SphericalBody(layers, 6371.0, name=:simple_earth)
```
"""
function SphericalBody(
    density,
    radius_km::Real;
    layer_boundaries::Union{Nothing, Vector{<:Real}}=nothing,
    name::Symbol=:unnamed
)
    T = Float64
    radius_natural = T(radius_km * units.km)

    # Handle different density specifications
    if density isa Vector || density isa Tuple
        # Layered density: vector of (density, boundary) tuples
        boundaries = T[0.0]
        density_funcs = DensityFunction{T, <:Any}[]

        for item in density
            d, b = item
            b_val = T(b)
            if !isempty(boundaries) && b_val <= boundaries[end]
                throw(ArgumentError("Layer boundaries must be strictly increasing. " *
                    "Got $b_val after $(boundaries[end])."))
            end
            push!(boundaries, b_val)
            df = DensityFunction(d, T(units.DENSITY_CONV))
            push!(density_funcs, df)
        end
    else
        # Single density value or function
        if isnothing(layer_boundaries)
            boundaries = T[0.0, 1.0]
        else
            boundaries = T.(layer_boundaries)
        end
        df = DensityFunction(density, T(units.DENSITY_CONV))
        density_funcs = [df]
    end

    # Compute average densities for each layer
    avg_densities = compute_average_densities(density_funcs, boundaries)

    SphericalBody{T}(radius_natural, boundaries, density_funcs, avg_densities, name)
end

"""
Compute average density in each layer via numerical integration.
"""
function compute_average_densities(
    density_funcs::Vector{<:DensityFunction{T}},
    boundaries::Vector{T}
) where T
    n_layers = length(density_funcs)
    avg_densities = Vector{T}(undef, n_layers)

    for i in 1:n_layers
        r_low = boundaries[i]
        r_high = boundaries[i + 1]

        if r_high > r_low
            integral, _ = quadgk(density_funcs[i], r_low, r_high)
            avg_densities[i] = T(integral / (r_high - r_low))
        else
            avg_densities[i] = density_funcs[i](r_low)
        end
    end

    return avg_densities
end

# Property accessors
"""Return radius in natural units (eV^-1)."""
radius(body::SphericalBody) = body.radius_natural

"""Return radius in kilometers."""
radius_km(body::SphericalBody) = body.radius_natural / units.km

"""Return length (same as radius for spherical body) in natural units."""
Base.length(body::SphericalBody) = body.radius_natural

"""Return layer boundaries."""
layer_boundaries(body::SphericalBody) = body.layer_boundaries

"""
    get_density(body::SphericalBody, r::Real; right::Bool=false)

Get density at normalized radius r ∈ [0, 1].

# Arguments
- `body`: The spherical body
- `r`: Normalized radius (0 = center, 1 = surface)
- `right`: If true, use right-side boundary for edge cases

# Returns
Density in natural units (eV^4)
"""
function get_density(body::SphericalBody{T}, r::Real; right::Bool=false) where T
    if r == 1
        layer_idx = length(body.density_functions)
    elseif r == 0
        layer_idx = 1
    else
        layer_idx = searchsortedlast(body.layer_boundaries, r)
        if right
            layer_idx = searchsortedfirst(body.layer_boundaries, r) - 1
        end
        layer_idx = clamp(layer_idx, 1, length(body.density_functions))
    end
    return body.density_functions[layer_idx](r)
end

"""
    get_average_density(body::SphericalBody, r::Real)

Get the average density of the layer containing radius r.

# Arguments
- `body`: The spherical body
- `r`: Normalized radius (0 = center, 1 = surface)

# Returns
Average density of the layer in natural units (eV^4)
"""
function get_average_density(body::SphericalBody{T}, r::Real) where T
    if !(0 <= r <= 1)
        throw(ArgumentError("r must be between 0 and 1, got $r"))
    end

    layer_idx = searchsortedlast(body.layer_boundaries, r, lt=<=)
    layer_idx = clamp(layer_idx, 1, length(body.average_densities))

    return body.average_densities[layer_idx]
end

function Base.show(io::IO, body::SphericalBody)
    n_layers = length(body.density_functions)
    r_km = radius_km(body)
    print(io, "SphericalBody(:$(body.name), radius=$(round(r_km, digits=1)) km, $n_layers layers)")
end
