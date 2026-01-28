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
            "density" => density_gcm3,
            "geometries" => [
                Dict(
                    "hierarchy" => i - 1,
                    "shape" => "sphere",
                    "origin" => [0.0, 0.0, 0.0],
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
                "hierarchy" => n_layers,
                "shape" => "sphere",
                "origin" => [0.0, 0.0, 0.0],
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
    generate_slab_config(body::AbstractSlabBody; cuts=default_cuts())

Generate a PROPOSAL JSON configuration for a slab body.

Uses box geometry to represent the slab.
"""
function generate_slab_config(body::AbstractSlabBody; cuts=default_cuts(), outer_buffer::Real=1e10)
    # Slab length in cm
    length_cm = body.length_natural / units.cm

    # Use a large transverse size (effectively infinite)
    transverse_size = 1e10  # cm

    sectors = []
    boundaries = layer_boundaries(body)

    z_offset = 0.0
    for i in 1:(length(boundaries) - 1)
        z_start = boundaries[i] * length_cm
        z_end = boundaries[i + 1] * length_cm
        layer_thickness = z_end - z_start

        # Get density and map to medium
        x_mid = (boundaries[i] + boundaries[i + 1]) / 2
        density = get_density(body, x_mid)
        medium_name = density_to_medium(density)

        sector = Dict(
            "medium" => medium_name,
            "geometries" => [
                Dict(
                    "hierarchy" => i - 1,
                    "shape" => "box",
                    "origin" => [0.0, 0.0, z_start + layer_thickness / 2],
                    "length" => transverse_size,
                    "width" => transverse_size,
                    "height" => layer_thickness
                )
            ]
        )

        push!(sectors, sector)
    end

    # Add air buffer beyond slab end so particles don't propagate
    # beyond the defined geometry (which causes PROPOSAL segfaults)
    n_layers = length(boundaries) - 1
    push!(sectors, Dict(
        "medium" => "Air",
        "geometries" => [
            Dict(
                "hierarchy" => n_layers,
                "shape" => "box",
                "origin" => [0.0, 0.0, length_cm + outer_buffer / 2],
                "length" => transverse_size,
                "width" => transverse_size,
                "height" => outer_buffer
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
