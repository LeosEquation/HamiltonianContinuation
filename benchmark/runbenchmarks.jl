using HamiltonianContinuation
using TaylorSeries
using BenchmarkTools

set_variables("x", numvars=5)

include("equilibria/continuation_1dof.jl")
println()
include("equilibria/continuation_2dof.jl")