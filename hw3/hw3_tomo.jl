# Problem 3 Tomography
# By: Koichi Kimoto

include(joinpath(@__DIR__, "..", "readclassjson.jl"))
using JSON
using LinearAlgebra
using Plots

data = readclassjson(joinpath(@__DIR__, "tomo_data.json"))

N = data["N"]
npixels = data["npixels"]
y = data["y"]
line_pixel_lengths = data["line_pixel_lengths"]

println("y length: $(length(y)), line_pixel_lengths size: $(size(line_pixel_lengths))")

x = line_pixel_lengths' \ y

X = reshape(x, npixels, npixels)

heatmap(X, yflip=true, aspect_ratio=:equal, color=:gist_gray,
        cbar=:none, framestyle=:none, title="Tomography Reconstruction")
savefig("tomography_reconstruction.png")