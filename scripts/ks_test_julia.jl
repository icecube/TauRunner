#!/usr/bin/env julia
# Run tau neutrinos through Earth at several angles and save outgoing energies
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "TauRunner.jl"))

using TauRunner
using DelimitedFiles

function main()
    n_events = parse(Int, ARGS[1])
    thetas_deg = parse.(Float64, split(ARGS[2], ","))
    energy = parse(Float64, ARGS[3])
    outdir = ARGS[4]
    seed = parse(Int, ARGS[5])

    earth = construct_earth()
    xs = CrossSections(CSMS)

    for θ_deg in thetas_deg
        θ_rad = deg2rad(θ_deg)
        energies = fill(energy, n_events)
        thetas = fill(θ_rad, n_events)

        results = run_mc(energies, thetas, earth, xs;
                         flavor=NuTau, seed=seed, secondaries=true, losses=true)

        # Collect final energies for particles that exited
        out_energies = Float64[]
        for r in results
            if r.position >= 1.0
                push!(out_energies, r.E_final)
            end
        end

        fname = joinpath(outdir, "julia_theta_$(θ_deg)_n_$(n_events).csv")
        writedlm(fname, out_energies, ',')
        println("θ=$(θ_deg)°: $(length(out_energies))/$(n_events) exited, saved to $fname")
    end
end

main()
