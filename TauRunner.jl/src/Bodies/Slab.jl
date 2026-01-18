"""
Slab body geometry for TauRunner.jl

Implements flat slab geometry for approximating particle propagation
through layered media without spherical geometry.
"""

"""
    LayeredConstantSlab{T} <: AbstractSlabBody

A slab body with constant density layers.

# Fields
- `length_natural::T`: Total length of the slab in natural units (eV^-1)
- `layer_boundaries::Vector{T}`: Normalized boundaries (0 to 1) for each layer
- `densities::Vector{T}`: Constant density for each layer in natural units (eV^4)
- `name::Symbol`: Name identifier for the body
"""
struct LayeredConstantSlab{T<:Real} <: AbstractSlabBody
    length_natural::T
    layer_boundaries::Vector{T}
    densities::Vector{T}
    name::Symbol
end

"""
    LayeredConstantSlab(layers, total_length_km; name=:slab)

Construct a LayeredConstantSlab from layer specifications.

# Arguments
- `layers`: Vector of (density_g_cm3, fractional_boundary) tuples
- `total_length_km`: Total length of the slab in kilometers
- `name`: Optional name identifier

# Examples
```julia
# Two-layer slab: rock (2.65 g/cm³) for first half, water (1.0 g/cm³) for second half
slab = LayeredConstantSlab([(2.65, 0.5), (1.0, 1.0)], 100.0)

# Single uniform layer
slab = LayeredConstantSlab([(2.65, 1.0)], 50.0, name=:rock_slab)
```
"""
function LayeredConstantSlab(
    layers::Vector{<:Tuple{Real, Real}},
    total_length_km::Real;
    name::Symbol=:slab
)
    T = Float64
    length_natural = T(total_length_km * units.km)

    boundaries = T[0.0]
    densities = T[]

    for (density, boundary) in layers
        push!(boundaries, T(boundary))
        push!(densities, T(density * units.DENSITY_CONV))
    end

    return LayeredConstantSlab{T}(length_natural, boundaries, densities, name)
end

"""
    LayeredConstantSlab(density, total_length_km; name=:uniform_slab)

Construct a uniform density slab.
"""
function LayeredConstantSlab(
    density::Real,
    total_length_km::Real;
    name::Symbol=:uniform_slab
)
    return LayeredConstantSlab([(density, 1.0)], total_length_km, name=name)
end

# Property accessors
Base.length(body::LayeredConstantSlab) = body.length_natural
layer_boundaries(body::LayeredConstantSlab) = body.layer_boundaries

"""
    get_density(body::LayeredConstantSlab, x::Real)

Get density at normalized position x ∈ [0, 1].
"""
function get_density(body::LayeredConstantSlab{T}, x::Real) where T
    if x == 1
        return body.densities[end]
    elseif x == 0
        return body.densities[1]
    end

    layer_idx = searchsortedlast(body.layer_boundaries, x)
    layer_idx = clamp(layer_idx, 1, length(body.densities))

    return body.densities[layer_idx]
end

"""
    get_average_density(body::LayeredConstantSlab, x::Real)

Get the average density of the layer containing position x.
For constant-density layers, this equals the local density.
"""
function get_average_density(body::LayeredConstantSlab{T}, x::Real) where T
    return get_density(body, x)
end

function Base.show(io::IO, body::LayeredConstantSlab)
    n_layers = length(body.densities)
    l_km = body.length_natural / units.km
    print(io, "LayeredConstantSlab(:$(body.name), length=$(round(l_km, digits=1)) km, $n_layers layers)")
end
