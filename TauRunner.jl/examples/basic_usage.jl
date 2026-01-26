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
# Example 7: Cross-Section Models (loading from JLD2)
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 7: Cross-Section Models (loading from JLD2)")
println("=" ^ 60)

println("\nTauRunner.jl supports two cross-section models:")
println("  - DIPOLE: Dipole model for neutrino-nucleon interactions")
println("  - CSMS:   Connolly-Sarkar-Mohanty-Sahu (perturbative QCD) model")

# Load both cross-section models from JLD2 files
println("\nLoading DIPOLE cross-sections from JLD2...")
xs_dipole = CrossSections(DIPOLE)
println("  Loaded: $xs_dipole")

println("\nLoading CSMS cross-sections from JLD2...")
xs_csms = CrossSections(CSMS)
println("  Loaded: $xs_csms")

# Conversion factor: natural units (eV^-2) to cm^2
# 1 eV^-1 = hbar*c / eV ≈ 1.97e-5 cm, so 1 eV^-2 ≈ 3.89e-10 cm^2
const EV2_TO_CM2 = (1.97326963e-5)^2  # (hbar*c in cm*eV)^2

# Test energies: 1 PeV, 10 PeV, 100 PeV, 1 EeV
test_energies_eV = [1e15, 1e16, 1e17, 1e18]  # in eV
test_energy_names = ["1 PeV", "10 PeV", "100 PeV", "1 EeV"]

println("\nTotal neutrino CC cross-sections (ν + N → l + X):")
println("-" ^ 55)
println("  Energy       |    DIPOLE         |    CSMS")
println("-" ^ 55)
for (E, name) in zip(test_energies_eV, test_energy_names)
    σ_dipole = total_cross_section(xs_dipole, E, :nu, :CC) * EV2_TO_CM2
    σ_csms = total_cross_section(xs_csms, E, :nu, :CC) * EV2_TO_CM2
    println("  $(rpad(name, 12)) |  $(round(σ_dipole / 1e-33, digits=2)) × 10⁻³³ cm²  |  $(round(σ_csms / 1e-33, digits=2)) × 10⁻³³ cm²")
end
println("-" ^ 55)

# Verify cross-sections are in expected range (physics sanity check)
# At 1 PeV, neutrino-nucleon CC cross-sections should be ~10^-33 cm^2
# At 1 EeV, should be ~10^-32 cm^2
println("\nPhysics sanity check:")
σ_1PeV = total_cross_section(xs_csms, 1e15, :nu, :CC) * EV2_TO_CM2
σ_1EeV = total_cross_section(xs_csms, 1e18, :nu, :CC) * EV2_TO_CM2
if 1e-35 < σ_1PeV < 1e-31 && 1e-34 < σ_1EeV < 1e-30
    println("  ✓ Cross-sections at 1 PeV and 1 EeV are in expected range")
else
    println("  ✗ WARNING: Cross-sections outside expected range!")
end

# Test differential cross-section
println("\nDifferential cross-section dσ/dz at 1 PeV (z = outgoing lepton energy fraction):")
E_test = 1e16  # 10 PeV
z_values = [0.1, 0.3, 0.5, 0.7, 0.9]
println("  z      |  dσ/dz (CSMS)")
for z in z_values
    dsigma = differential_cross_section(xs_csms, E_test, z, :nu, :CC) * EV2_TO_CM2
    println("  $(z)    |  $(round(dsigma / 1e-33, digits=3)) × 10⁻³³ cm²")
end

# =============================================================================
# Example 8: Monte Carlo Neutrino Propagation
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 8: Monte Carlo Neutrino Propagation")
println("=" ^ 60)

using Random

# Set up the simulation components
earth = construct_earth()
xs = CrossSections(CSMS)

println("\nSingle particle propagation:")
println("-" ^ 40)

# Create a tau neutrino with 100 PeV energy
E_initial = 1e17  # 100 PeV in eV
theta = π/4       # 45° nadir angle

# Create the track and propagator
track = Chord(theta)
clp = SphericalBodyPropagator(earth)

# Create a particle
# Particle(type, energy, position, cross_sections; secondaries=true, losses=true)
particle = Particle(NuTau, E_initial, 0.0, xs; secondaries=true, losses=true)

println("  Initial: $(particle.id) at E = $(E_initial/GeV) GeV")
println("  Track: Chord at θ = $(round(theta * 180/π, digits=1))°")

# Propagate with a fixed random seed for reproducibility
rng = MersenneTwister(42)
propagate!(particle, track, earth, clp; rng=rng)

println("\n  Final state:")
println("    Particle: $(particle.id)")
println("    Energy: $(round(particle.energy/GeV, digits=1)) GeV")
println("    Position: $(round(particle.position, digits=4)) (1.0 = exited)")
println("    Survived: $(particle.survived)")
println("    CC interactions: $(particle.n_cc)")
println("    NC interactions: $(particle.n_nc)")
println("    Decays: $(particle.n_decay)")

# =============================================================================
# Example 9: Batch Monte Carlo with run_mc()
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 9: Batch Monte Carlo Simulation")
println("=" ^ 60)

# Run multiple events at once using run_mc()
n_events = 10000
E_batch = 1e18  # 100 PeV
theta_batch = 0.0  # Nadir (straight through Earth center)

# Create arrays of energies and angles
energies = fill(E_batch, n_events)
thetas = fill(theta_batch, n_events)

println("\nRunning $n_events events at E = $(E_batch/1e15) PeV, θ = 0° (nadir)...")
println("(This may take a moment on first run while PROPOSAL builds tables)")

# Run the Monte Carlo
results = run_mc(energies, thetas, earth, xs;
                 flavor=NuTau,
                 seed=12345,
                 losses=true,
                 secondaries=true)

# Analyze results
n_exited = count(r -> r.position >= 1.0, results)
n_nutau = count(r -> r.id == Int(NuTau), results)
n_tau = count(r -> r.id == Int(Tau), results)

# Energy statistics for exiting particles
exit_energies = [r.E_final for r in results if r.position >= 1.0]
mean_E_out = isempty(exit_energies) ? 0.0 : sum(exit_energies) / length(exit_energies)

# Interaction statistics
total_cc = sum(r.n_cc for r in results)
total_nc = sum(r.n_nc for r in results)

println("\nResults:")
println("  Events exited: $n_exited / $n_events ($(round(100*n_exited/n_events, digits=1))%)")
println("  Final particles: ν_τ=$n_nutau, τ=$n_tau")
println("  Mean exit energy: $(round(mean_E_out/GeV, digits=1)) GeV")
println("  Mean CC interactions: $(round(total_cc/n_events, digits=2))")
println("  Mean NC interactions: $(round(total_nc/n_events, digits=2))")

# =============================================================================
# Example 10: Energy-dependent Transmission
# =============================================================================

println("\n" * "=" ^ 60)
println("Example 10: Energy-dependent Transmission")
println("=" ^ 60)

println("\nComputing transmission probability vs energy at θ = 60°...")

n_per_energy = 50
theta_test = π/3  # 60°
test_Es = [1e16, 1e17, 1e18]  # 10 PeV, 100 PeV, 1 EeV

println("\n  Energy     | Exit rate | Mean E_out/E_in")
println("  " * "-" ^ 45)

for E in test_Es
    local energies = fill(E, n_per_energy)
    local thetas = fill(theta_test, n_per_energy)

    local results = run_mc(energies, thetas, earth, xs;
                           flavor=NuTau, seed=999, losses=true)

    n_exit = count(r -> r.position >= 1.0, results)
    exit_frac = n_exit / n_per_energy

    # Mean energy ratio for exiting particles
    ratios = [r.E_final / r.E_initial for r in results if r.position >= 1.0]
    mean_ratio = isempty(ratios) ? 0.0 : sum(ratios) / length(ratios)

    E_name = E >= 1e18 ? "$(E/1e18) EeV" : "$(E/1e15) PeV"
    println("  $(rpad(E_name, 10)) |   $(round(100*exit_frac, digits=1))%   |   $(round(mean_ratio, digits=3))")
end

println("\n" * "=" ^ 60)
println("Examples complete!")
println("=" ^ 60)
