"""
Earth model for TauRunner.jl

Implements the Preliminary Reference Earth Model (PREM) with 9 layers.
"""

# PREM polynomial coefficients for each layer
# Density = a0 + a1*r + a2*r^2 + a3*r^3  (in g/cm^3, r normalized)
const PREM_PARAMS = [
    (13.0885,   0.0,    -8.8381,   0.0),      # Inner core
    (12.5815,  -1.2638, -3.6426,  -5.5281),   # Outer core
    (7.9565,   -6.4761,  5.5283,  -3.0807),   # Lower mantle
    (5.3197,   -1.4836,  0.0,      0.0),      # Transition zone 1
    (11.2494,  -8.0298,  0.0,      0.0),      # Transition zone 2
    (7.1089,   -3.8045,  0.0,      0.0),      # Transition zone 3
    (2.6910,    0.6924,  0.0,      0.0),      # LVZ/LID
    (2.9,       0.0,     0.0,      0.0),      # Crust (oceanic)
    (2.6,       0.0,     0.0,      0.0),      # Crust (continental)
]

# Layer boundaries in km from center (standard PREM)
const PREM_BOUNDARIES_KM = [1221.5, 3480, 5701, 5771, 5971, 6151, 6346.6, 6356, 6368]
const PREM_RADIUS_KM = 6368.0

"""
    prem_density(r, params)

Evaluate PREM polynomial density at normalized radius r.
"""
function prem_density(r::Real, params::NTuple{4, <:Real})
    a0, a1, a2, a3 = params
    return a0 + a1*r + a2*r^2 + a3*r^3
end

"""
Create a density function closure for a given set of PREM parameters.
"""
function make_prem_density_func(params::NTuple{4, <:Real})
    return r -> prem_density(r, params)
end

"""
    construct_earth(; layers::Vector{Tuple{Real, Real}}=Tuple{Real, Real}[])

Construct a PREM Earth model.

# Arguments
- `layers`: Optional list of (thickness_km, density_g_cm3) tuples to add as
  constant density layers on top of the PREM model. For example, a 4 km
  ocean layer would be `[(4.0, 1.0)]`.

# Returns
A `SphericalBody` representing the PREM Earth, optionally with additional layers.

# Examples
```julia
# Standard PREM Earth
earth = construct_earth()

# PREM Earth with 4 km ocean
earth_ocean = construct_earth(layers=[(4.0, 1.0)])

# PREM Earth with 3 km ice and 4 km ocean
earth_ice = construct_earth(layers=[(3.0, 0.92), (4.0, 1.0)])
```
"""
function construct_earth(; layers::Vector{<:Tuple{Real, Real}}=Tuple{Real, Real}[])
    r_tot = PREM_RADIUS_KM
    boundary_list_km = collect(Float64, PREM_BOUNDARIES_KM)
    params_list = collect(PREM_PARAMS)

    # Add additional surface layers
    for (thickness, density) in layers
        r_tot += thickness
        push!(params_list, (density, 0.0, 0.0, 0.0))
        push!(boundary_list_km, r_tot)
    end

    # Normalize boundaries to [0, 1]
    normalized_boundaries = boundary_list_km ./ r_tot

    # Create density function for each layer
    earth_layers = Tuple{Function, Float64}[]
    for (params, boundary) in zip(params_list, normalized_boundaries)
        density_func = make_prem_density_func(params)
        push!(earth_layers, (density_func, boundary))
    end

    return SphericalBody(earth_layers, r_tot, name=:PREM_earth)
end

"""
    construct_uniform_earth(density_g_cm3::Real=5.5)

Construct a uniform density Earth for testing/comparison.

# Arguments
- `density_g_cm3`: Uniform density in g/cm^3 (default: 5.5, roughly Earth average)
"""
function construct_uniform_earth(density_g_cm3::Real=5.5)
    return SphericalBody(density_g_cm3, PREM_RADIUS_KM, name=:uniform_earth)
end
