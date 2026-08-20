module HamiltonianContinuation

using TaylorIntegration
using TaylorSeries
using LinearAlgebra
using PeriodicSchurDecompositions
using Printf

# Integration
include("integration/integration.jl")

# Objects
include("objects/objects.jl")

# Equilibria
include("equilibria/equilibria.jl")

# Limit Points
include("limitpoint/limitpoint.jl")

# Branch Points
include("branchpoint/branchpoint.jl")

# Orbits
include("orbits/orbits.jl")

# Flip Points
include("flippoint/flippoint.jl")


export BPContinuation,
    EquilibriaContinuation,
    FPContinuation,
    LPContinuation,
    OrbitContinuation

export HamiltonianBranchPoint,
    HamiltonianBranchPointBranch,
    HamiltonianLimitPoint,
    HamiltonianLimitPointBranch,
    HamiltonianEquilibriumBranch,
    HamiltonianFlipPoint,
    HamiltonianFlipPointBranch,
    HamiltonianPeriodicBranch

end