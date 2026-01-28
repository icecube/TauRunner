"""
Cross-section model definitions and data loading.
"""

"""
Cross-section model types.
"""
@enum XSModel begin
    DIPOLE  # Dipole model
    CSMS    # Perturbative QCD model (Cooper-Sarkar, Mertsch, Sarkar)
end

"""
    CrossSections{T} <: AbstractCrossSections

Container for neutrino cross-section splines.

# Fields
- `model::XSModel`: The cross-section model (DIPOLE or CSMS)
- `total_splines::Dict`: Splines for total cross-sections
- `differential_splines::Dict`: Splines for differential cross-sections
"""
struct CrossSections{T<:Real} <: AbstractCrossSections
    model::XSModel
    # Keys: (nutype, target, interaction) where:
    #   nutype: :nu or :nubar
    #   target: :p (proton) or :n (neutron)
    #   interaction: :CC or :NC
    total_splines::Dict{Tuple{Symbol, Symbol, Symbol}, Any}
    differential_splines::Dict{Tuple{Symbol, Symbol, Symbol}, Any}
end

"""
    CrossSections(model::XSModel; datapath::String=default_data_path())

Load cross-section data for the specified model.

# Arguments
- `model`: Cross-section model (DIPOLE or CSMS)
- `datapath`: Path to data files (default: package data directory)

# Examples
```julia
xs = CrossSections(CSMS)
xs = CrossSections(DIPOLE)
```
"""
function CrossSections(model::XSModel; datapath::Union{String, Nothing}=nothing)
    T = Float64

    if isnothing(datapath)
        datapath = default_data_path()
    end

    total_splines = Dict{Tuple{Symbol, Symbol, Symbol}, Any}()
    differential_splines = Dict{Tuple{Symbol, Symbol, Symbol}, Any}()

    model_name = model == DIPOLE ? "dipole" : "csms"
    model_path = joinpath(datapath, "cross_sections", model_name)

    # Load splines for each combination
    for nutype in (:nu, :nubar)
        for target in (:p, :n)
            for interaction in (:CC, :NC)
                key = (nutype, target, interaction)

                # Load total cross-section spline
                total_file = joinpath(model_path, "$(nutype)_$(target)_sigma_$(interaction).jld2")
                if isfile(total_file)
                    total_splines[key] = load_spline(total_file)
                else
                    error("Required cross-section file not found: $total_file")
                end

                # Load differential cross-section spline (2D: energy x z)
                diff_file = joinpath(model_path, "$(nutype)_$(target)_dsigma_$(interaction).jld2")
                if isfile(diff_file)
                    differential_splines[key] = load_spline_2d(diff_file)
                end
            end
        end
    end

    return CrossSections{T}(model, total_splines, differential_splines)
end

"""
Default path to cross-section data files.
"""
function default_data_path()
    # Look for data directory relative to package
    pkg_dir = dirname(dirname(@__DIR__))
    data_path = joinpath(pkg_dir, "data")

    if !isdir(data_path)
        @warn "Data directory not found at $data_path. Cross-sections will need to be loaded manually."
    end

    return data_path
end

"""
    total_cross_section(xs::CrossSections, energy, nutype, interaction; proton_fraction=0.5)

Compute the total neutrino cross-section at given energy.

# Arguments
- `xs`: CrossSections object
- `energy`: Neutrino energy in eV
- `nutype`: :nu or :nubar (or "nu"/"nubar" string)
- `interaction`: :CC or :NC (or "CC"/"NC" string)
- `proton_fraction`: Fraction of target that is protons (default: 0.5 for isoscalar)

# Returns
Total cross-section in natural units (eV^-2)
"""
function total_cross_section(
    xs::CrossSections{T},
    energy::Real,
    nutype::Union{Symbol, String},
    interaction::Union{Symbol, String};
    proton_fraction::Real=0.5
) where T
    nutype_sym = nutype isa String ? Symbol(nutype) : nutype
    interaction_sym = interaction isa String ? Symbol(interaction) : interaction

    neutron_fraction = 1.0 - proton_fraction

    key_p = (nutype_sym, :p, interaction_sym)
    key_n = (nutype_sym, :n, interaction_sym)

    log_E = log(energy)

    # Get cross-sections from splines (stored as log(sigma) vs log(E))
    if haskey(xs.total_splines, key_p) && haskey(xs.total_splines, key_n)
        sigma_p = exp(xs.total_splines[key_p](log_E))
        sigma_n = exp(xs.total_splines[key_n](log_E))
        return T(proton_fraction * sigma_p + neutron_fraction * sigma_n)
    else
        error("Cross-section data not available for $key_p or $key_n")
    end
end

"""
    differential_cross_section(xs::CrossSections, energy, z, nutype, interaction; proton_fraction=0.5)

Compute the differential cross-section dσ/dz at given energy and inelasticity.

# Arguments
- `xs`: CrossSections object
- `energy`: Neutrino energy in eV
- `z`: Energy fraction retained by outgoing lepton (0 < z < 1)
- `nutype`: :nu or :nubar
- `interaction`: :CC or :NC
- `proton_fraction`: Fraction of target that is protons

# Returns
Differential cross-section dσ/dz in natural units
"""
function differential_cross_section(
    xs::CrossSections{T},
    energy::Real,
    z::Union{Real, AbstractVector{<:Real}},
    nutype::Union{Symbol, String},
    interaction::Union{Symbol, String};
    proton_fraction::Real=0.5
) where T
    nutype_sym = nutype isa String ? Symbol(nutype) : nutype
    interaction_sym = interaction isa String ? Symbol(interaction) : interaction

    neutron_fraction = 1.0 - proton_fraction

    key_p = (nutype_sym, :p, interaction_sym)
    key_n = (nutype_sym, :n, interaction_sym)

    log_E = log(energy)

    if haskey(xs.differential_splines, key_p) && haskey(xs.differential_splines, key_n)
        # Differential splines are 2D: f(log_E, z), stored as log(dsigma/dE)
        # To get dσ/dz, divide by energy (matching Python implementation)
        if z isa AbstractVector
            dsigma_p = [exp(xs.differential_splines[key_p](log_E, zi)) / energy for zi in z]
            dsigma_n = [exp(xs.differential_splines[key_n](log_E, zi)) / energy for zi in z]
            return T.(proton_fraction .* dsigma_p .+ neutron_fraction .* dsigma_n)
        else
            dsigma_p = exp(xs.differential_splines[key_p](log_E, z)) / energy
            dsigma_n = exp(xs.differential_splines[key_n](log_E, z)) / energy
            return T(proton_fraction * dsigma_p + neutron_fraction * dsigma_n)
        end
    else
        error("Differential cross-section data not available for $key_p or $key_n")
    end
end

function Base.show(io::IO, xs::CrossSections)
    print(io, "CrossSections($(xs.model), $(length(xs.total_splines)) total, $(length(xs.differential_splines)) differential)")
end
