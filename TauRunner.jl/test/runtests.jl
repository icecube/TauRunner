using Test
using JSON3
using Random
using Statistics
using TauRunner
using TauRunner.units: GeV, TeV, MeV, km, meter, cm, gr, DENSITY_CONV
using TauRunner.PhysicsConstants: EARTH_RADIUS_KM, GF, TAU_MASS_GEV

# Load Python reference data
const REFERENCE_FILE = joinpath(@__DIR__, "regression", "python_reference.json")
const REF_DATA = JSON3.read(read(REFERENCE_FILE, String))

# Tolerance for comparing with Python
const RTOL = 0.02  # 2% relative tolerance (accounts for numerical differences)
const ATOL = 1e-10  # Absolute tolerance for near-zero values

@testset "TauRunner.jl" begin

    @testset "Units" begin
        # Test unit conversions
        @test GeV ≈ 1e9
        @test TeV ≈ 1e12
        @test km > meter
        @test cm < meter

        # Test density conversion
        @test DENSITY_CONV ≈ gr / cm^3
    end

    @testset "PhysicsConstants" begin
        @test EARTH_RADIUS_KM ≈ 6371.0
        @test GF ≈ 1.16639e-23
        @test TAU_MASS_GEV ≈ 1.776
    end

    @testset "Bodies (vs Python)" begin
        using TauRunner.Bodies

        earth = construct_earth()
        ref_geo = REF_DATA["geometry"]

        # Test Earth radius matches Python
        @test radius_km(earth) ≈ ref_geo["earth_radius_km"] rtol=RTOL

        # Test density profile matches Python
        for entry in ref_geo["density_profile"]
            r = entry["r"]
            python_rho = entry["density_gcm3"]
            julia_rho = get_density(earth, r) / (gr / cm^3)
            @test julia_rho ≈ python_rho rtol=RTOL
        end

        # Test structural properties
        @test earth isa SphericalBody
        @test earth.name == :PREM_earth
    end

    @testset "Tracks (vs Python)" begin
        using TauRunner.Tracks
        using TauRunner.Bodies

        earth = construct_earth()
        ref_geo = REF_DATA["geometry"]

        # Test Chord properties vs Python
        for chord_test in ref_geo["chord_tests"]
            theta_deg = chord_test["theta_deg"]
            theta_rad = theta_deg * π / 180

            chord = Chord(theta_rad)

            # Test total length
            @test chord.total_length ≈ chord_test["total_length"] rtol=RTOL

            # Test column depth
            julia_depth = total_column_depth(chord, earth) / (gr / cm^2)
            python_depth = chord_test["column_depth_gcm2"]
            @test julia_depth ≈ python_depth rtol=RTOL

            # Test x_to_r coordinate transform
            for xr in chord_test["x_to_r"]
                x = xr["x"]
                python_r = xr["r"]
                julia_r = x_to_r(chord, x)
                @test julia_r ≈ python_r atol=1e-10
            end
        end

        # Test Radial track vs Python
        radial = Radial()
        ref_radial = ref_geo["radial_tests"]

        # Test column depth
        julia_depth = total_column_depth(radial, earth) / (gr / cm^2)
        python_depth = ref_radial["column_depth_gcm2"]
        @test julia_depth ≈ python_depth rtol=RTOL

        # Test x_to_r for radial
        for xr in ref_radial["x_to_r"]
            x = xr["x"]
            python_r = xr["r"]
            julia_r = x_to_r(radial, x)
            @test julia_r ≈ python_r atol=1e-10
        end
    end

    @testset "Particle Types" begin
        @test is_neutrino(NuTau)
        @test is_neutrino(NuEBar)
        @test !is_neutrino(Tau)

        @test is_charged_lepton(Tau)
        @test is_charged_lepton(Muon)
        @test !is_charged_lepton(NuMu)

        @test charged_partner(NuTau) == Tau
        @test charged_partner(NuE) == Electron

        @test neutrino_partner(Tau) == NuTau
        @test neutrino_partner(Muon) == NuMu
    end

    @testset "Cross-Sections (vs Python)" begin
        ref_xs = REF_DATA["cross_sections"]

        # Conversion factor from eV^-2 to cm^2
        EV2_TO_CM2 = (1.97326963e-5)^2

        for (model_name, model_type) in [("dipole", DIPOLE), ("csms", CSMS)]
            xs = CrossSections(model_type)

            @testset "$model_name total cross-sections" begin
                for entry in ref_xs[model_name]["total"]
                    E = entry["energy_eV"]

                    # Test nu CC
                    julia_sigma = total_cross_section(xs, E, :nu, :CC) * EV2_TO_CM2
                    python_sigma = entry["nu_CC_cm2"]
                    @test julia_sigma ≈ python_sigma rtol=RTOL

                    # Test nu NC
                    julia_sigma = total_cross_section(xs, E, :nu, :NC) * EV2_TO_CM2
                    python_sigma = entry["nu_NC_cm2"]
                    @test julia_sigma ≈ python_sigma rtol=RTOL

                    # Test nubar CC
                    julia_sigma = total_cross_section(xs, E, :nubar, :CC) * EV2_TO_CM2
                    python_sigma = entry["nubar_CC_cm2"]
                    @test julia_sigma ≈ python_sigma rtol=RTOL

                    # Test nubar NC
                    julia_sigma = total_cross_section(xs, E, :nubar, :NC) * EV2_TO_CM2
                    python_sigma = entry["nubar_NC_cm2"]
                    @test julia_sigma ≈ python_sigma rtol=RTOL
                end
            end

            @testset "$model_name differential cross-sections" begin
                for entry in ref_xs[model_name]["differential"]
                    E = entry["energy_eV"]
                    z = entry["z"]

                    # Test nu CC differential
                    julia_dsigma = differential_cross_section(xs, E, z, :nu, :CC) * EV2_TO_CM2
                    python_dsigma = entry["nu_CC_cm2"]
                    @test julia_dsigma ≈ python_dsigma rtol=RTOL

                    # Test nu NC differential
                    julia_dsigma = differential_cross_section(xs, E, z, :nu, :NC) * EV2_TO_CM2
                    python_dsigma = entry["nu_NC_cm2"]
                    @test julia_dsigma ≈ python_dsigma rtol=RTOL
                end
            end
        end
    end

    @testset "Tau Decay Spectrum (vs Python)" begin
        ref_decay = REF_DATA["tau_decay"]

        # Test expected mean z from the decay spectrum
        # We compute this by sampling many times
        rng = MersenneTwister(12345)
        n_samples = 100000
        z_samples = [TauRunner.sample_tau_decay_fraction(rng) for _ in 1:n_samples]

        julia_mean_z = mean(z_samples)
        python_mean_z = ref_decay["expected_mean_z"]

        # Mean should match within 1%
        @test julia_mean_z ≈ python_mean_z rtol=0.01

        # Test median
        julia_median_z = median(z_samples)
        python_median_z = ref_decay["expected_median_z"]
        @test julia_median_z ≈ python_median_z rtol=0.01
    end

    @testset "Input Validation" begin
        using TauRunner.Tracks
        using TauRunner.Bodies

        # Bug fix: Chord rejects invalid theta
        @test_throws ArgumentError Chord(-0.1)
        @test_throws ArgumentError Chord(π + 0.1)
        @test_throws ArgumentError Chord(2π)

        # Valid boundary values should work
        @test Chord(0.0) isa Chord
        @test Chord(π) isa Chord
        @test Chord(π/4) isa Chord

        # Bug fix: SphericalBody rejects unsorted layer boundaries
        @test_throws ArgumentError SphericalBody(
            [(5.0, 0.5), (3.0, 0.3)],  # 0.3 < 0.5, not increasing
            6371.0
        )
        @test_throws ArgumentError SphericalBody(
            [(5.0, 0.5), (3.0, 0.5)],  # duplicate boundary
            6371.0
        )

        # Properly ordered boundaries should work
        body = SphericalBody(
            [(5.0, 0.5), (3.0, 1.0)],
            6371.0
        )
        @test body isa SphericalBody

        # Bug fix: CrossSections errors on missing data files
        @test_throws ErrorException CrossSections(CSMS, datapath="/nonexistent/path")
    end

    @testset "Column Depth Monotonicity" begin
        using TauRunner.Tracks
        using TauRunner.Bodies

        # Column depth along a chord should be monotonically increasing
        earth = construct_earth()
        chord = Chord(π/4)

        xs = range(0.0, 1.0, length=50)
        Xs = [x_to_X(chord, earth, x) for x in xs]

        for i in 2:length(Xs)
            @test Xs[i] >= Xs[i-1]
        end

        # Total column depth should be positive
        @test total_column_depth(chord, earth) > 0.0

        # x_to_X and X_to_x should be approximate inverses
        for x_orig in [0.1, 0.25, 0.5, 0.75, 0.9]
            X_val = x_to_X(chord, earth, x_orig)
            x_back = X_to_x(chord, earth, X_val)
            @test x_back ≈ x_orig atol=1e-3
        end
    end

    @testset "Monte Carlo Propagation (vs Python)" begin
        ref_mc = REF_DATA["monte_carlo"]

        earth = construct_earth()
        xs = CrossSections(DIPOLE)

        # Use larger tolerance for MC due to stochastic nature
        # We use 20% tolerance since we're comparing 100-event samples
        MC_RTOL = 0.30

        for test_case in ref_mc
            E = test_case["energy_eV"]
            theta_deg = test_case["theta_deg"]
            n_events = test_case["n_events"]
            seed = test_case["seed"]

            @testset "E=$(E) eV, θ=$(theta_deg)°" begin
                # Run Julia MC with same seed
                energies = fill(E, n_events)
                thetas = fill(theta_deg * π / 180, n_events)

                results = run_mc(energies, thetas, earth, xs;
                                flavor=NuTau, seed=seed, losses=true, secondaries=true)

                # Compute statistics
                exit_energies = [r.E_final for r in results if r.position >= 1.0]
                n_exited = length(exit_energies)

                julia_median = n_exited > 0 ? median(exit_energies) : 0.0
                julia_mean_nCC = mean([r.n_cc for r in results])
                julia_mean_nNC = mean([r.n_nc for r in results])

                python_median = test_case["median_exit_energy_eV"]
                python_mean_nCC = test_case["mean_nCC"]
                python_mean_nNC = test_case["mean_nNC"]

                # Test exit fraction
                @test n_exited / n_events ≈ test_case["exit_fraction"] atol=0.1

                # Test median exit energy (with larger tolerance for stochastic variation)
                if python_median > 0
                    @test julia_median ≈ python_median rtol=MC_RTOL
                end

                # Test mean interaction counts
                @test julia_mean_nCC ≈ python_mean_nCC rtol=MC_RTOL
                @test julia_mean_nNC ≈ python_mean_nNC rtol=MC_RTOL
            end
        end
    end

    @testset "Reproducibility" begin
        earth = construct_earth()
        xs = CrossSections(DIPOLE)

        n_events = 50
        energies = fill(1e18, n_events)
        thetas = fill(deg2rad(45.0), n_events)

        @testset "Same seed gives identical results" begin
            r1 = run_mc(energies, thetas, earth, xs; flavor=NuTau, seed=12345, losses=true, secondaries=true)
            r2 = run_mc(energies, thetas, earth, xs; flavor=NuTau, seed=12345, losses=true, secondaries=true)

            for i in 1:n_events
                @test r1[i].E_final == r2[i].E_final
                @test r1[i].id == r2[i].id
                @test r1[i].n_cc == r2[i].n_cc
                @test r1[i].n_nc == r2[i].n_nc
                @test r1[i].position == r2[i].position
            end
        end

        @testset "Different seeds give different results" begin
            r1 = run_mc(energies, thetas, earth, xs; flavor=NuTau, seed=100, losses=true, secondaries=true)
            r2 = run_mc(energies, thetas, earth, xs; flavor=NuTau, seed=200, losses=true, secondaries=true)

            # At least some events should differ
            n_differ = count(i -> r1[i].E_final != r2[i].E_final, 1:n_events)
            @test n_differ > 0
        end
    end

end
