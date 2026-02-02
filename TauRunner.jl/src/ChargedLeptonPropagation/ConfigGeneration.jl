"""
Configuration generation for PROPOSAL.jl propagators.

Generates JSON configuration files from TauRunner body definitions.
PROPOSAL uses JSON configs to define geometry, media, and energy cuts.
"""

# PROPOSAL units: position in cm, energy in MeV

# Default path for PROPOSAL interpolation tables
# This is set to a subdirectory within the TauRunner.jl package
const DEFAULT_PROPOSAL_TABLES_PATH = abspath(joinpath(@__DIR__, "..", "..", "data", "proposal_tables"))

"""
    get_proposal_tables_path()

Get the path for PROPOSAL interpolation tables.
Returns the default path within TauRunner.jl/data/proposal_tables,
or the value of PROPOSAL_TABLES_PATH environment variable if set.
"""
function get_proposal_tables_path()
    path = get(ENV, "PROPOSAL_TABLES_PATH", DEFAULT_PROPOSAL_TABLES_PATH)
    # Ensure the directory exists
    mkpath(path)
    return abspath(path)
end

"""
Map TauRunner density (in natural units eV^4) to PROPOSAL medium name.
This is approximate - we select the closest standard medium.
"""
function density_to_medium(density_natural::Real)
    # Convert from natural units to g/cm³
    density_gcm3 = density_natural / units.DENSITY_CONV

    # Standard medium densities (g/cm³) from PROPOSAL
    media = [
        ("Air", 0.001225),
        ("Ice", 0.917),
        ("Water", 1.0),
        ("StandardRock", 2.65),
        ("FrejusRock", 2.74),
        ("Iron", 7.874),
        ("Lead", 11.35),
    ]

    # Find closest match
    best_medium = "StandardRock"
    best_diff = Inf

    for (name, rho) in media
        diff = abs(rho - density_gcm3)
        if diff < best_diff
            best_diff = diff
            best_medium = name
        end
    end

    return best_medium
end

"""
Default energy cut settings for PROPOSAL.
"""
function default_cuts()
    return Dict(
        "e_cut" => 1e14,        # Effectively infinite energy cut (matches Python TauRunner's np.inf)
        "v_cut" => 1e-3,        # Relative energy cut (matches Python TauRunner)
        "cont_rand" => true     # Continuous randomization
    )
end

"""
    generate_sphere_config(body::AbstractSphericalBody; cuts=default_cuts(), outer_buffer=1e10)

Generate a PROPOSAL JSON configuration for a spherical body.

The body is represented as concentric spherical shells, each with the
medium that best matches the layer's density. An outer "Air" layer is added
to handle particles at/near the surface.

# Arguments
- `body`: A SphericalBody from TauRunner
- `cuts`: Energy cut settings (optional)
- `outer_buffer`: Additional radius beyond body surface for air layer (cm, default: 1e10)

# Returns
A Dict that can be serialized to JSON for PROPOSAL.
"""
function generate_sphere_config(body::AbstractSphericalBody; cuts=default_cuts(), outer_buffer::Real=1e10)
    # Get body radius in cm (PROPOSAL units)
    radius_cm = radius(body) / units.cm

    # Build sectors for each layer
    sectors = []
    boundaries = layer_boundaries(body)

    # Use sphere origin (0,0,1) to match Python TauRunner convention.
    # The 1 cm z-offset breaks exact symmetry at shell boundaries, preventing
    # PROPOSAL's IsInside check from hitting the degenerate case where a
    # particle lies exactly on a boundary (determinant == 0 is excluded).
    sphere_origin = [0.0, 0.0, 1.0]

    for i in 1:(length(boundaries) - 1)
        r_inner = boundaries[i] * radius_cm
        r_outer = boundaries[i + 1] * radius_cm

        # Get average density in this layer (g/cm³)
        r_mid = (boundaries[i] + boundaries[i + 1]) / 2
        density = get_density(body, r_mid)
        density_gcm3 = density / units.DENSITY_CONV

        # Use StandardRock for all layers with explicit density override
        # (matches Python TauRunner which uses density_homogeneous)
        sector = Dict(
            "medium" => "StandardRock",
            "density_distribution" => Dict("type" => "homogeneous", "mass_density" => density_gcm3),
            "geometries" => [
                Dict(
                    "hierarchy" => 0,
                    "shape" => "sphere",
                    "origin" => sphere_origin,
                    "outer_radius" => r_outer,
                    "inner_radius" => r_inner
                )
            ]
        )

        push!(sectors, sector)
    end

    # Add outer air buffer layer so particles near the surface don't
    # propagate beyond the defined geometry (which causes PROPOSAL segfaults)
    n_layers = length(boundaries) - 1
    push!(sectors, Dict(
        "medium" => "Air",
        "geometries" => [
            Dict(
                "hierarchy" => 0,
                "shape" => "sphere",
                "origin" => sphere_origin,
                "outer_radius" => radius_cm + outer_buffer,
                "inner_radius" => radius_cm
            )
        ]
    ))

    config = Dict(
        "global" => Dict(
            "cuts" => cuts,
            "tables_path" => get_proposal_tables_path()
        ),
        "sectors" => sectors
    )

    return config
end

"""
    generate_uniform_sphere_config(radius_cm::Real, medium::String; cuts=default_cuts())

Generate a simple uniform sphere configuration.

# Arguments
- `radius_cm`: Radius in centimeters
- `medium`: PROPOSAL medium name (e.g., "StandardRock", "Ice", "Water")
- `cuts`: Energy cut settings
"""
function generate_uniform_sphere_config(radius_cm::Real, medium::String; cuts=default_cuts())
    config = Dict(
        "global" => Dict(
            "cuts" => cuts,
            "tables_path" => get_proposal_tables_path()
        ),
        "sectors" => [
            Dict(
                "medium" => medium,
                "geometries" => [
                    Dict(
                        "hierarchy" => 0,
                        "shape" => "sphere",
                        "origin" => [0.0, 0.0, 0.0],
                        "outer_radius" => Float64(radius_cm)
                    )
                ]
            )
        ]
    )
    return config
end

"""
    generate_slab_layer_config(density_gcm3::Real, medium_name::String; cuts=default_cuts())

Generate a PROPOSAL JSON configuration for a single slab layer.

Uses a large enclosing sphere (radius 1e20 cm) centered at the origin,
matching the Python TauRunner approach. The particle is always placed at
the origin and propagated for the layer thickness. This avoids all
boundary-related PROPOSAL segfaults ("No sector defined at particle
position") that occur with fitted box geometries when a particle lands
exactly on a sector boundary due to floating-point arithmetic.
"""
function generate_slab_layer_config(density_gcm3::Real, medium_name::String; cuts=default_cuts())
    config = Dict(
        "global" => Dict(
            "cuts" => cuts,
            "tables_path" => get_proposal_tables_path()
        ),
        "sectors" => [
            Dict(
                "medium" => medium_name,
                "density_distribution" => Dict("type" => "homogeneous", "mass_density" => density_gcm3),
                "geometries" => [
                    Dict(
                        "hierarchy" => 0,
                        "shape" => "sphere",
                        "origin" => [0.0, 0.0, 0.0],
                        "outer_radius" => 1e20
                    )
                ]
            )
        ]
    )

    return config
end

"""
    write_config(config::Dict, filepath::String)

Write a PROPOSAL configuration to a JSON file.
"""
function write_config(config::Dict, filepath::String)
    open(filepath, "w") do io
        JSON3.pretty(io, config)
    end
    return filepath
end

"""
    create_temp_config(config::Dict)

Write config to a temporary file and return the path.
"""
function create_temp_config(config::Dict)
    filepath = tempname() * ".json"
    write_config(config, filepath)
    return filepath
end
