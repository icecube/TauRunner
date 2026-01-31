#!/usr/bin/env julia
# Compare Julia and Python TauRunner outputs using KS test (primary particles only)
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "TauRunner.jl"))

using DelimitedFiles, HypothesisTests, Statistics, Printf

n_events = 10000
thetas = [5.0, 15.0, 30.0, 45.0, 60.0, 80.0]
outdir = joinpath(@__DIR__, "ks_output")

println("="^80)
println("KS Test: Julia vs Python TauRunner (primary particles only)")
println("Energy = 1e17 eV, N = $n_events per angle")
println("="^80)
@printf("%-8s  %8s  %8s  %12s  %10s  %8s\n",
        "θ [deg]", "N_jl", "N_py", "KS stat", "p-value", "Result")
println("-"^80)

for θ in thetas
    jl_file = joinpath(outdir, "julia_theta_$(θ)_n_$(n_events).csv")
    py_file = joinpath(outdir, "python_theta_$(θ)_n_$(n_events).csv")

    jl_data = readdlm(jl_file, Float64) |> vec
    py_data = readdlm(py_file, Float64) |> vec

    # Filter: exited particles have E > 0 (non-exited recorded as 0)
    jl_exited = filter(x -> x > 0, jl_data)
    py_exited = filter(x -> x > 0, py_data)

    ks = ApproximateTwoSampleKSTest(jl_exited, py_exited)
    pval = pvalue(ks)
    result = pval < 0.05 ? "REJECT" : "PASS"

    @printf("%-8.1f  %8d  %8d  %12.6f  %10.4f  %8s\n",
            θ, length(jl_exited), length(py_exited), ks.δ, pval, result)
end
println("="^80)
