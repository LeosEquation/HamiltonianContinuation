using Test
using HamiltonianContinuation
using TaylorSeries

set_variables("x", numvars=6)

include("equilibria/continuation_1dof.jl")
include("equilibria/continuation_2dof.jl")