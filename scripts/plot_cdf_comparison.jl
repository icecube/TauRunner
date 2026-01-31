using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "TauRunner.jl"))

using DelimitedFiles, HypothesisTests, Statistics
using Plots
gr()

E0 = 1e17
thetas = [5.0, 15.0, 30.0, 45.0, 60.0, 80.0, 89.0]
outdir = joinpath(@__DIR__, "ks_output")

plts = []

for θ in thetas
    jl_file = joinpath(outdir, "julia_theta_$(θ)_n_10000.csv")
    py_file = joinpath(outdir, "python_theta_$(θ)_n_10000.csv")
    if !isfile(py_file) || !isfile(jl_file)
        println("Skipping θ=$(θ)° — missing data")
        continue
    end

    jl = readdlm(jl_file, Float64) |> vec
    py = readdlm(py_file, Float64) |> vec

    jl_pos = filter(x -> x > 0, jl)
    py_pos = filter(x -> x > 0, py)

    jl_log = sort(log10.(jl_pos))
    py_log = sort(log10.(py_pos))

    jl_cdf = (1:length(jl_log)) ./ length(jl_log)
    py_cdf = (1:length(py_log)) ./ length(py_log)

    ks = ApproximateTwoSampleKSTest(jl_pos, py_pos)
    pval = round(pvalue(ks), digits=3)

    sp = plot(jl_log, jl_cdf, label="Julia", color=:blue, lw=1.5,
              xlabel="log₁₀(E / eV)", ylabel="CDF",
              title="θ = $(θ)° (p = $(pval))", legend=:bottomright)
    plot!(sp, py_log, py_cdf, label="Python", color=:red, lw=1.5, ls=:dash)
    push!(plts, sp)
end

fig = plot(plts..., layout=(length(plts), 1), size=(600, 300 * length(plts)), dpi=150,
           left_margin=5Plots.mm, bottom_margin=3Plots.mm)
savefig(fig, joinpath(outdir, "cdf_comparison.png"))
println("Saved to $(joinpath(outdir, "cdf_comparison.png"))")
