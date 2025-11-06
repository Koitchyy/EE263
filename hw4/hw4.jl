# Problem 4 Optimal Trade Off
# By: Koichi Kimoto

include(joinpath(@__DIR__, "..", "readclassjson.jl"))
using JSON
using LinearAlgebra
using Plots

data = readclassjson("curve_smoothing.json")
n = data["n"]
f = data["f"]

A = zeros(n - 2, n)
for i in range(1, n - 2)
    A[i, i] = 1
    A[i, i + 1] = -2
    A[i, i + 2] = 1
end

mu_vals = [0, 0.0001, 0.001, 0.1, 10^6]

for mu in mu_vals
    g = inv((1/n)*I + mu * (n^4/(n-2)) * A' * A) * (1/n) * f
    plot([f g],
    label=["f" "g with μ=$mu"],
    linestyle=[:dot :solid],
    linewidth=2)    
    savefig("f_vs_g_µ=$mu.png")
end

mu_sweep = logrange(10^-4, 10^4, 1000)

c = Float64[]
d = Float64[]

for mu in mu_sweep
    # find minimal g
    g = inv((1/n)*I + mu * (n^4/(n-2)) * A' * A) * (1/n) * f
    # get coord of (c, d)
    push!(c, (n^4/(n - 2)) * norm(A * g)^2)
    push!(d, (1/n) * norm(f - g)^2)
end

plot(d, c, title="Optimal Trade-off Curve",
label="Pareto Curve",
xlabel = "d",
ylabel = "c",
linewidth=2)

savefig("tradeoff.png")





