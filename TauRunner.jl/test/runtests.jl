using Test
using TauRunner
using TauRunner.units: GeV, TeV, MeV, km, meter, cm, gr, DENSITY_CONV
using TauRunner.PhysicsConstants: EARTH_RADIUS_KM, GF, TAU_MASS_GEV

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

    @testset "Bodies" begin
        using TauRunner.Bodies

        # Test Earth construction
        earth = construct_earth()
        @test earth isa SphericalBody
        @test earth.name == :PREM_earth

        # Test radius
        @test radius_km(earth) ≈ 6368.0 atol=1.0

        # Test density at different radii
        # Core should be denser than surface
        core_density = get_density(earth, 0.0)
        surface_density = get_density(earth, 1.0)
        @test core_density > surface_density

        # Test average density
        avg_core = get_average_density(earth, 0.1)
        @test avg_core > 0
    end

    @testset "Tracks" begin
        using TauRunner.Tracks
        using TauRunner.Bodies

        # Test Chord construction
        chord = Chord(0.0)  # Through center
        @test chord.theta ≈ 0.0
        @test chord.total_length ≈ 2.0  # Diameter

        # Test coordinate transforms
        @test x_to_r(chord, 0.0) ≈ 1.0  # At surface
        @test x_to_r(chord, 0.5) ≈ 0.0 atol=1e-10  # At center
        @test x_to_r(chord, 1.0) ≈ 1.0  # Back at surface

        # Test Radial track
        radial = Radial()
        @test x_to_r(radial, 0.0) ≈ 1.0  # At surface
        @test x_to_r(radial, 1.0) ≈ 0.0  # At center

        # Test column depth calculation
        earth = construct_earth()
        depth = total_column_depth(chord, earth)
        @test depth > 0
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

end
