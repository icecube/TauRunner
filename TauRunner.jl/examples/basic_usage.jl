#=
TauRunner.jl Basic Usage Examples

This file demonstrates basic usage of the TauRunner.jl package for
propagating ultra-high-energy neutrinos through matter.

NOTE: Full Monte Carlo propagation requires cross-section data to be loaded.
These examples focus on geometry and column depth calculations which work
without the full cross-section infrastructure.
=#

using TauRunner
using TauRunner.units: GeV, TeV, km, cm, gr

# =============================================================================
# Example 1: Create Earth with PREM density model
# =============================================================================

println("=" ^ 60)
println("Example 1: Earth Model")
println("=" ^ 60)

# Construct Earth with the PREM (Preliminary Reference Earth Model)
earth = construct_earth()

println("Earth model: ", earth.name)
println("Radius: ", radius_km(earth), " km")

# Check density at different radii (r is normalized: 0=center, 1=surface)
println("\nDensity profile:")
for r in [0.0, 0.2, 0.5, 0.8, 1.0]
    ρ = get_density(earth, r) / (gr / cm^3)
    println("  r = $r R_Earth: ρ = $(round(ρ, digits=2)) g/cm³")
end

# =============================================================================
# Example 2: Chord trajectories through Earth
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 2: Chord Trajectories")
println("=" ^ 60)

# A chord is defined by the nadir angle θ (in radians)
# θ = 0: straight through the center (diameter)
# θ = π/2: grazing the surface

println("\nChord properties at different nadir angles:")
for θ_deg in [0.0, 30.0, 60.0, 85.0, 89.0]
    θ = θ_deg * π / 180
    chord = Chord(θ)

    # Total path length (in units of Earth radius)
    L = chord.total_length
    L_km = L * radius_km(earth)

    println("  θ = $(θ_deg)°: path length = $(round(L_km, digits=1)) km")
end

# =============================================================================
# Example 3: Column depth calculations
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 3: Column Depth")
println("=" ^ 60)

# Column depth is the integrated density along the path: ∫ρ·dl
# Important for determining interaction probability

println("\nColumn depth for different trajectories:")
for θ_deg in [0.0, 30.0, 60.0, 85.0, 89.0]
    θ = θ_deg * π / 180
    chord = Chord(θ)

    # Calculate total column depth
    X = total_column_depth(chord, earth)
    X_gcm2 = X / (gr / cm^2)

    println("  θ = $(θ_deg)°: X = $(round(X_gcm2 / 1e9, digits=2)) × 10⁹ g/cm²")
end

# =============================================================================
# Example 4: Coordinate transformations along a track
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 4: Track Coordinates")
println("=" ^ 60)

# Create a chord through Earth at 45° nadir angle
chord = Chord(45.0 * π / 180)

println("\nPosition along chord (θ = 45°):")
println("  x: fractional position along track (0=entry, 1=exit)")
println("  r: normalized radius (0=center, 1=surface)")

for x in [0.0, 0.25, 0.5, 0.75, 1.0]
    r = x_to_r(chord, x)
    println("  x = $x → r = $(round(r, digits=3))")
end

# =============================================================================
# Example 5: Radial trajectory
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 5: Radial Trajectory")
println("=" ^ 60)

# Radial trajectory goes from surface to center
radial = Radial()

println("\nRadial trajectory (surface to center):")
for x in [0.0, 0.25, 0.5, 0.75, 1.0]
    r = x_to_r(radial, x)
    println("  x = $x → r = $(round(r, digits=3))")
end

# Column depth for radial trajectory
X_radial = total_column_depth(radial, earth)
println("\nRadial column depth: $(round(X_radial / (gr / cm^2) / 1e9, digits=2)) × 10⁹ g/cm²")

# =============================================================================
# Example 6: Particle types
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 6: Particle Types")
println("=" ^ 60)

println("\nNeutrino types:")
for p in [NuE, NuMu, NuTau, NuEBar, NuMuBar, NuTauBar]
    partner = charged_partner(p)
    println("  $p (PDG=$(Int(p))) → charged partner: $partner")
end

println("\nCharged lepton types:")
for p in [Electron, Muon, Tau]
    partner = neutrino_partner(p)
    println("  $p (PDG=$(Int(p))) → neutrino partner: $partner")
end

# =============================================================================
# Example 7: Using different cross-section models
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 7: Cross-Section Models (requires data files)")
println("=" ^ 60)

println("""
TauRunner.jl supports two cross-section models:
  - DIPOLE: Dipole model for neutrino-nucleon interactions
  - CSMS:   Connolly-Sarkar-Mohanty-Sahu model

To use cross-sections, load the data files:

    using TauRunner
    xs = CrossSections(DIPOLE)  # or CSMS

    # Total cross-section at energy E
    σ = total_cross_section(xs, E, :nu, :CC)

    # Differential cross-section
    dσ_dy = differential_cross_section(xs, E_in, y, :nu, :CC)
""")

println("\n" * "=" ^ 60)
println("Examples complete!")
println("=" ^ 60)
