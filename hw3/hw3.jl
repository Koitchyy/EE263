# Problem 3 Sensor Integrity
# By: Koichi Kimoto

include(joinpath(@__DIR__, "..", "readclassjson.jl"))
using JSON
using LinearAlgebra

A = [1  2  1 
     1 -1 -2
    -2  1  3
     1 -1 -2
     1  1  0]

println(rank(A))

N = nullspace(A')
B = N'
println("B = ")
display(B)

# Check for By = 0
x = [-1 0 1]'
y = A * x
println("Should be a 3-vector with entries ~ 0")
display(B * y)
