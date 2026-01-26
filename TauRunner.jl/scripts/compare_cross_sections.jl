#!/usr/bin/env julia
"""
Comprehensive comparison of Julia TauRunner.jl with Python TauRunner.

This script outputs values for all aspects covered by the regression tests:
- Units and conversion factors
- Physics constants
- Bodies/Geometry (Earth radius, density profile)
- Tracks (Chord properties, x_to_r, column depth)
- Particle types
- Cross-sections (total and differential)

Run from the TauRunner.jl directory:
    julia --project=. scripts/compare_cross_sections.jl
"""

using TauRunner
using TauRunner.units: GeV, TeV, MeV, km, meter, cm, gr, sec, DENSITY_CONV
using TauRunner.PhysicsConstants: EARTH_RADIUS_KM, GF, TAU_MASS_GEV,
    MUON_MASS_GEV, PROTON_MASS_GEV, NEUTRON_MASS_GEV, TAU_LIFETIME_SEC
using Printf
using Random

# Conversion factor: natural units (eV^-2) to cm^2
const EV2_TO_CM2 = (1.97326963e-5)^2

function print_section(title)
    println("\n" * "=" ^ 70)
    println(title)
    println("=" ^ 70)
end

function print_subsection(title)
    println("\n--- $title ---")
end

# =============================================================================
# UNITS COMPARISON
# =============================================================================
function compare_units()
    print_section("UNITS AND CONVERSION FACTORS")

    println("\nEnergy units (eV per unit):")
    println("  GeV  = $(GeV)")
    println("  TeV  = $(TeV)")
    println("  MeV  = $(MeV)")

    println("\nDistance units (eV^-1 per unit):")
    println("  km    = $(km)")
    println("  meter = $(meter)")
    println("  cm    = $(cm)")

    println("\nMass units (eV per unit):")
    println("  gr = $(gr)")

    println("\nTime units (eV^-1 per unit):")
    println("  sec = $(sec)")

    println("\nDerived quantities:")
    println("  DENSITY_CONV (gr/cm^3) = $(DENSITY_CONV)")
end

# =============================================================================
# PHYSICS CONSTANTS COMPARISON
# =============================================================================
function compare_constants()
    print_section("PHYSICS CONSTANTS")

    println("\nGeometry:")
    println("  EARTH_RADIUS_KM = $(EARTH_RADIUS_KM) km")

    println("\nFundamental constants:")
    println("  GF (Fermi constant) = $(GF) eV^-2")

    println("\nParticle masses (GeV):")
    println("  TAU_MASS_GEV      = $(TAU_MASS_GEV)")
    println("  MUON_MASS_GEV     = $(MUON_MASS_GEV)")
    println("  PROTON_MASS_GEV   = $(PROTON_MASS_GEV)")
    println("  NEUTRON_MASS_GEV  = $(NEUTRON_MASS_GEV)")

    println("\nLifetimes:")
    println("  TAU_LIFETIME_SEC = $(TAU_LIFETIME_SEC) s")
end

# =============================================================================
# BODIES/GEOMETRY COMPARISON
# =============================================================================
function compare_geometry()
    print_section("BODIES / GEOMETRY")

    earth = construct_earth()
    println("\nEarth model: $(earth.name)")
    println("  Radius (km): $(radius_km(earth))")
    println("  Radius (eV^-1): $(radius(earth))")

    print_subsection("Density Profile")
    println("  r (fractional) | density (g/cm³) | density (eV^4)")
    println("  " * "-" ^ 50)

    r_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for r in r_values
        rho_ev4 = get_density(earth, r)
        rho_gcm3 = rho_ev4 / DENSITY_CONV
        @printf("  %.1f              %8.4f          %.6e\n", r, rho_gcm3, rho_ev4)
    end

    print_subsection("Average Density")
    println("  r (fractional) | avg density (g/cm³)")
    println("  " * "-" ^ 40)

    for r in [0.1, 0.3, 0.5, 0.7, 0.9]
        avg_rho = get_average_density(earth, r) / DENSITY_CONV
        @printf("  %.1f              %8.4f\n", r, avg_rho)
    end
end

# =============================================================================
# TRACKS COMPARISON
# =============================================================================
function compare_tracks()
    print_section("TRACKS")

    earth = construct_earth()

    # Test angles (degrees)
    theta_tests = [0.0, 30.0, 45.0, 60.0, 90.0]

    print_subsection("Chord Properties")

    for theta_deg in theta_tests
        theta_rad = theta_deg * π / 180
        chord = Chord(theta_rad)

        println("\n  θ = $(theta_deg)°:")
        println("    Total length (dimensionless): $(chord.total_length)")

        # Column depth
        col_depth = total_column_depth(chord, earth)
        col_depth_gcm2 = col_depth / (gr / cm^2)
        @printf("    Column depth: %.6e g/cm²\n", col_depth_gcm2)
        @printf("    Column depth: %.6e eV³\n", col_depth)

        # x_to_r mapping
        println("    x_to_r mapping:")
        x_values = [0.0, 0.25, 0.5, 0.75, 1.0]
        for x in x_values
            r = x_to_r(chord, x)
            @printf("      x=%.2f -> r=%.6f\n", x, r)
        end
    end

    print_subsection("Radial Track")
    radial = Radial()
    println("  x_to_r mapping (radial):")
    for x in [0.0, 0.25, 0.5, 0.75, 1.0]
        r = x_to_r(radial, x)
        @printf("    x=%.2f -> r=%.6f\n", x, r)
    end
end

# =============================================================================
# PARTICLE TYPES COMPARISON
# =============================================================================
function compare_particle_types()
    print_section("PARTICLE TYPES")

    print_subsection("Particle IDs (PDG)")
    particles = [NuE, NuMu, NuTau, NuEBar, NuMuBar, NuTauBar,
                 Electron, Muon, Tau, Positron, AntiMuon, AntiTau]

    println("  Particle     | PDG ID | is_neutrino | is_charged_lepton")
    println("  " * "-" ^ 60)

    for p in particles
        @printf("  %-12s | %6d | %-11s | %s\n",
                p, Int(p), is_neutrino(p), is_charged_lepton(p))
    end

    print_subsection("Charged Partners (neutrino -> charged lepton)")
    for nu in [NuE, NuMu, NuTau, NuEBar, NuMuBar, NuTauBar]
        partner = charged_partner(nu)
        println("  $(nu) -> $(partner)")
    end

    print_subsection("Neutrino Partners (charged lepton -> neutrino)")
    for cl in [Electron, Muon, Tau, Positron, AntiMuon, AntiTau]
        partner = neutrino_partner(cl)
        println("  $(cl) -> $(partner)")
    end
end

# =============================================================================
# CROSS-SECTIONS COMPARISON
# =============================================================================
function compare_cross_sections()
    print_section("CROSS-SECTIONS")

    # Load cross-sections
    println("\nLoading CSMS cross-sections...")
    xs_csms = CrossSections(CSMS)
    println("  Loaded: $(xs_csms)")

    println("\nLoading DIPOLE cross-sections...")
    xs_dipole = CrossSections(DIPOLE)
    println("  Loaded: $(xs_dipole)")

    # Test energies
    test_energies_eV = [1e14, 1e15, 1e16, 1e17, 1e18, 1e19]
    test_energy_names = ["100 TeV", "1 PeV", "10 PeV", "100 PeV", "1 EeV", "10 EeV"]

    # -------------------------------------------------------------------------
    # Total Cross-Sections
    # -------------------------------------------------------------------------
    print_subsection("Total Cross-Sections (CSMS)")

    for (interaction, int_name) in [(:CC, "CC"), (:NC, "NC")]
        println("\n  $int_name interactions:")
        println("  Energy       | nu (cm²)             | nubar (cm²)")
        println("  " * "-" ^ 55)

        for (E, name) in zip(test_energies_eV, test_energy_names)
            sigma_nu = total_cross_section(xs_csms, E, :nu, interaction) * EV2_TO_CM2
            sigma_nubar = total_cross_section(xs_csms, E, :nubar, interaction) * EV2_TO_CM2
            @printf("  %-12s | %.4e        | %.4e\n", name, sigma_nu, sigma_nubar)
        end
    end

    print_subsection("Total Cross-Sections (DIPOLE)")

    for (interaction, int_name) in [(:CC, "CC"), (:NC, "NC")]
        println("\n  $int_name interactions:")
        println("  Energy       | nu (cm²)             | nubar (cm²)")
        println("  " * "-" ^ 55)

        for (E, name) in zip(test_energies_eV, test_energy_names)
            sigma_nu = total_cross_section(xs_dipole, E, :nu, interaction) * EV2_TO_CM2
            sigma_nubar = total_cross_section(xs_dipole, E, :nubar, interaction) * EV2_TO_CM2
            @printf("  %-12s | %.4e        | %.4e\n", name, sigma_nu, sigma_nubar)
        end
    end

    # -------------------------------------------------------------------------
    # Differential Cross-Sections
    # -------------------------------------------------------------------------
    print_subsection("Differential Cross-Sections (CSMS)")

    E_test = 1e16  # 10 PeV
    z_values = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]

    println("\n  E = 1e16 eV (10 PeV)")

    for (interaction, int_name) in [(:CC, "CC"), (:NC, "NC")]
        println("\n  $int_name dσ/dz:")
        println("  z       | nu (cm²)             | nubar (cm²)")
        println("  " * "-" ^ 50)

        for z in z_values
            dsigma_nu = differential_cross_section(xs_csms, E_test, z, :nu, interaction) * EV2_TO_CM2
            dsigma_nubar = differential_cross_section(xs_csms, E_test, z, :nubar, interaction) * EV2_TO_CM2
            @printf("  %.2f    | %.4e        | %.4e\n", z, dsigma_nu, dsigma_nubar)
        end
    end

    print_subsection("Differential Cross-Sections (DIPOLE)")

    println("\n  E = 1e16 eV (10 PeV)")

    for (interaction, int_name) in [(:CC, "CC"), (:NC, "NC")]
        println("\n  $int_name dσ/dz:")
        println("  z       | nu (cm²)             | nubar (cm²)")
        println("  " * "-" ^ 50)

        for z in z_values
            dsigma_nu = differential_cross_section(xs_dipole, E_test, z, :nu, interaction) * EV2_TO_CM2
            dsigma_nubar = differential_cross_section(xs_dipole, E_test, z, :nubar, interaction) * EV2_TO_CM2
            @printf("  %.2f    | %.4e        | %.4e\n", z, dsigma_nu, dsigma_nubar)
        end
    end
end

# =============================================================================
# MONTE CARLO COMPARISON
# =============================================================================
function compare_monte_carlo()
    print_section("MONTE CARLO PROPAGATION")

    # Load cross-sections and create Earth
    xs = CrossSections(CSMS)
    earth = construct_earth()

    # -------------------------------------------------------------------------
    # Interaction Depths
    # -------------------------------------------------------------------------
    print_subsection("Interaction Depths (mean column depth to interaction)")

    test_energies = [1e15, 1e16, 1e17, 1e18]
    energy_names = ["1 PeV", "10 PeV", "100 PeV", "1 EeV"]

    println("\n  CC interaction depth (g/cm²):")
    println("  Energy       | nu                  | nubar")
    println("  " * "-" ^ 50)

    for (E, name) in zip(test_energies, energy_names)
        # Create particles to get interaction depths
        p_nu = Particle(NuTau, E, 0.0, xs)
        p_nubar = Particle(NuTauBar, E, 0.0, xs)

        depth_nu = get_interaction_depth(p_nu, :CC) / (gr / cm^2)
        depth_nubar = get_interaction_depth(p_nubar, :CC) / (gr / cm^2)

        @printf("  %-12s | %.4e        | %.4e\n", name, depth_nu, depth_nubar)
    end

    println("\n  NC interaction depth (g/cm²):")
    println("  Energy       | nu                  | nubar")
    println("  " * "-" ^ 50)

    for (E, name) in zip(test_energies, energy_names)
        p_nu = Particle(NuTau, E, 0.0, xs)
        p_nubar = Particle(NuTauBar, E, 0.0, xs)

        depth_nu = get_interaction_depth(p_nu, :NC) / (gr / cm^2)
        depth_nubar = get_interaction_depth(p_nubar, :NC) / (gr / cm^2)

        @printf("  %-12s | %.4e        | %.4e\n", name, depth_nu, depth_nubar)
    end

    # -------------------------------------------------------------------------
    # Single Particle Propagation (deterministic with fixed seed)
    # -------------------------------------------------------------------------
    print_subsection("Single Particle Propagation (seed=42)")

    test_cases = [
        (1e17, 0.0, "100 PeV, nadir"),
        (1e17, π/4, "100 PeV, 45°"),
        (1e18, 0.0, "1 EeV, nadir"),
        (1e18, π/3, "1 EeV, 60°"),
    ]

    println("\n  Initial conditions -> Final state")
    println("  " * "-" ^ 70)

    for (E, theta, desc) in test_cases
        # Use fixed seed for reproducibility
        rng = Random.MersenneTwister(42)

        track = Chord(theta)
        clp = SphericalBodyPropagator(earth)
        particle = Particle(NuTau, E, 0.0, xs; secondaries=true, losses=true)

        propagate!(particle, track, earth, clp; rng=rng)

        E_out_GeV = particle.energy / GeV
        @printf("  %s:\n", desc)
        @printf("    Final particle: %s\n", particle.id)
        @printf("    E_out: %.4e GeV\n", E_out_GeV)
        @printf("    Position: %.6f\n", particle.position)
        @printf("    Survived: %s\n", particle.survived)
        @printf("    n_CC: %d, n_NC: %d, n_decay: %d\n",
                particle.n_cc, particle.n_nc, particle.n_decay)
        println()
    end

    # -------------------------------------------------------------------------
    # Statistical MC Run
    # -------------------------------------------------------------------------
    print_subsection("Statistical Monte Carlo (N=1000, seed=12345)")

    n_events = 1000
    seed = 12345

    mc_tests = [
        (1e17, 0.0, "100 PeV, nadir"),
        (1e17, π/4, "100 PeV, 45°"),
        (1e18, 0.0, "1 EeV, nadir"),
        (1e18, π/4, "1 EeV, 45°"),
    ]

    for (E, theta, desc) in mc_tests
        println("\n  $desc (N=$n_events):")

        # Create arrays for run_mc
        energies = fill(E, n_events)
        thetas = fill(theta, n_events)

        results = run_mc(energies, thetas, earth, xs;
                        flavor=NuTau, seed=seed, losses=true, secondaries=true)

        # Compute statistics
        n_survived = count(r -> r.position >= 1.0 && r.E_final > 0, results)
        n_nutau = count(r -> r.id == Int(NuTau), results)
        n_tau = count(r -> r.id == Int(Tau), results)
        n_other = n_events - n_nutau - n_tau

        # Energy statistics for survivors
        survivor_energies = [r.E_final for r in results if r.position >= 1.0 && r.E_final > 0]
        nutau_energies = [r.E_final for r in results if r.id == Int(NuTau) && r.position >= 1.0]
        tau_energies = [r.E_final for r in results if r.id == Int(Tau) && r.position >= 1.0]

        # CC/NC statistics
        total_cc = sum(r.n_cc for r in results)
        total_nc = sum(r.n_nc for r in results)
        mean_cc = total_cc / n_events
        mean_nc = total_nc / n_events

        @printf("    Exited: %d / %d (%.1f%%)\n", n_survived, n_events, 100*n_survived/n_events)
        @printf("    Final particles: ν_τ=%d, τ=%d, other=%d\n", n_nutau, n_tau, n_other)
        @printf("    Mean n_CC: %.3f, Mean n_NC: %.3f\n", mean_cc, mean_nc)

        if length(nutau_energies) > 0
            mean_E_nutau = sum(nutau_energies) / length(nutau_energies)
            @printf("    ν_τ <E_out>: %.4e eV (N=%d)\n", mean_E_nutau, length(nutau_energies))
        end

        if length(tau_energies) > 0
            mean_E_tau = sum(tau_energies) / length(tau_energies)
            @printf("    τ <E_out>: %.4e eV (N=%d)\n", mean_E_tau, length(tau_energies))
        end
    end

    # -------------------------------------------------------------------------
    # Tau Regeneration Test
    # -------------------------------------------------------------------------
    print_subsection("Tau Regeneration Chain (seed=42)")

    println("\n  Testing ν_τ -> τ -> ν_τ regeneration:")

    # Run multiple events to find one with regeneration
    rng = Random.MersenneTwister(42)
    n_regen = 0
    n_test = 100

    for i in 1:n_test
        track = Chord(0.0)  # Through center
        clp = SphericalBodyPropagator(earth)
        particle = Particle(NuTau, 1e18, 0.0, xs)
        propagate!(particle, track, earth, clp; rng=rng)

        if particle.n_cc > 0 && particle.n_decay > 0
            n_regen += 1
        end
    end

    @printf("    Events with regeneration (n_CC>0 && n_decay>0): %d / %d\n", n_regen, n_test)

    # Show one detailed example
    rng = Random.MersenneTwister(42)
    track = Chord(0.0)
    clp = SphericalBodyPropagator(earth)
    particle = Particle(NuTau, 1e18, 0.0, xs)
    propagate!(particle, track, earth, clp; rng=rng)

    println("\n    Example event (seed=42, E=1 EeV, nadir):")
    @printf("      n_CC=%d, n_NC=%d, n_decay=%d\n",
            particle.n_cc, particle.n_nc, particle.n_decay)
    @printf("      Final: %s at E=%.4e eV\n", particle.id, particle.energy)
    @printf("      Position: %.6f, Survived: %s\n", particle.position, particle.survived)
end

# =============================================================================
# MAIN
# =============================================================================
function main()
    println("=" ^ 70)
    println("TauRunner.jl Comprehensive Comparison Output")
    println("=" ^ 70)
    println("\nThis output is designed to be compared with Python TauRunner.")
    println("Run the corresponding Python script to generate comparison values.")

    compare_units()
    compare_constants()
    compare_geometry()
    compare_tracks()
    compare_particle_types()
    compare_cross_sections()
    compare_monte_carlo()

    print_section("COMPARISON COMPLETE")
    println("\nCompare these values with the Python output to verify consistency.")
end

main()
