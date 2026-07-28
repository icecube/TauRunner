using Test
using TauRunner
using TauRunner.units: GeV, eV, km
using TauRunner.Bodies: LayeredConstantSlab
using TauRunner.Tracks: SlabTrack
using TauRunner.Particles: Particle
using TauRunner.ChargedLeptonPropagationModule: SlabPropagator,
    propagate_charged_lepton!, seed_proposal!

@testset "Decay in flight is detected" begin
    # A 1e16 eV tau in standard rock has a boosted decay length of ~0.5 km, so
    # over a 10 km slab it decays with probability ~1 -- and it decays in flight,
    # with an energy millions of times its rest mass. The old rest-mass check
    # (E <= m_tau) never fires for such a decay, so this exercises the has_decay path.
    slab = LayeredConstantSlab(2.65, 10.0)
    prop = SlabPropagator(slab)
    track = SlabTrack()
    xs = CrossSections(DIPOLE)

    E0 = 1e16  # eV

    seed_proposal!(42)
    n_decayed = 0
    for trial in 1:5
        particle = Particle(Tau, E0, 0.0, xs)
        propagate_charged_lepton!(prop, particle, track)

        if particle.n_decay > 0
            n_decayed += 1
            # Decay happened strictly inside the slab
            @test 0.0 < particle.decay_position < 1.0
            # The tau became a tau neutrino carrying a sizeable fraction of E0,
            # not a ~rest-mass-energy remnant
            @test particle.id == NuTau
            @test particle.energy > 1e-3 * E0
        end
    end

    # Decay is essentially certain over 20 decay lengths; require all 5
    @test n_decayed == 5
end
