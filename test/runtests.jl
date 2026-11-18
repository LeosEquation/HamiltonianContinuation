using Test
using HamiltonianContinuation
using TaylorSeries

variables!("x", numvars=6, order=2)
# jetspace = JetSpace(
#     order=2,  # usa aquí el orden que necesites
#     variables=["x₁", "x₂", "x₃", "x₄", "x₅", "x₆"]
# )

include("equilibria/continuation_1dof.jl")
include("equilibria/continuation_2dof.jl")