using HamiltonianContinuation
using TaylorSeries
using BenchmarkTools

variables!("x", numvars=5, order=2)

include("equilibria/continuation_1dof.jl")
println()
include("equilibria/continuation_2dof.jl")