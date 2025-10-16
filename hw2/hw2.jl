# Problem 2
# By: Koichi Kimoto
# Checks if the test light response is in the span of the phosphor color response

include(joinpath(@__DIR__, "..", "readclassjson.jl"))
using JSON
using LinearAlgebra

########################### Part (c) ###################################

data = readclassjson(joinpath(@__DIR__, "color_perception_data.json"))
wavelength = data["wavelength"]
r = data["R_phosphor"]
g = data["G_phosphor"]
b = data["B_phosphor"]
l = data["L_coefficients"]
m = data["M_coefficients"]
s = data["S_coefficients"]
test_light = data["test_light"]

A = [l m s]'

c_test = A * test_light
c_r = A * r
c_g = A * g
c_b = A * b

B = hcat(c_r, c_g, c_b)  # columns are spanning vectors               
x = B \ c_test
in_span = isapprox(B*x, c_test; rtol=0, atol=1e-10)
println("Is it in the span?: $in_span")
if in_span
    println("Weights (in order of R, G, B): $x")
else
    println("Error: $(norm(B*x - c_test))")
end

########################### Part (d) ##########################################
# Import source light spectra
sunlight_I = data["sunlight"]
tungsten_I = data["tungsten"]

function nullspace(M; tol=nothing)
    F = svd(M; full=true)
    s = F.S
    if tol === nothing
        tol = maximum(size(M)) * eps(eltype(s)) * maximum(s)
    end
    r = count(>(tol), s) 
    return F.Vt'[:, r+1:end]
end

# Build matrices mapping reflectance -> LMS under each illuminant
M_sun = A * Diagonal(sunlight_I)
M_tun = A * Diagonal(tungsten_I)

# Nullspace for tungsten (metamers under tungsten)
N_tun = nullspace(M_tun)
direction = N_tun[:, 1] 

baseline = fill(0.5, length(direction))

# Find maximum alpha such that r0 +- alpha*d is in [0, 1]^n
function sym_alpha(r0, d)
    alphas_plus = Float64[]
    alphas_minus = Float64[]
    for i in eachindex(d)
        if d[i] > 0
            push!(alphas_plus, (1.0 - r0[i]) / d[i])
            push!(alphas_minus, r0[i] / d[i])
        elseif d[i] < 0
            push!(alphas_plus, r0[i] / (-d[i]))
            push!(alphas_minus, (1.0 - r0[i]) / (-d[i]))
        end
    end
    return 0.9 * min(minimum(alphas_plus), minimum(alphas_minus))
end

alpha = sym_alpha(baseline, direction)
objA_reflect = baseline + alpha * direction
objB_reflect = baseline - alpha * direction

# Verify: identical under tungsten, different under sunlight
c_t_objA = M_tun * objA_reflect
c_t_objB = M_tun * objB_reflect
c_s_objA = M_sun * objA_reflect
c_s_objB = M_sun * objB_reflect

println("||c_t_objA - c_t_objB|| = $(norm(c_t_objA - c_t_objB))  (≈ 0)")
println("||c_s_objA - c_s_objB|| = $(norm(c_s_objA - c_s_objB))  (> 0)")
println("Example reflectances (first 5 entries):")
println("  objA: ", round.(objA_reflect[1:5]; digits=4))
println("  objB: ", round.(objB_reflect[1:5]; digits=4))





