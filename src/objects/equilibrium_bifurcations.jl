# Complete hierarchy for bifurcation points in general and Hamiltonian systems

# General abstract type for all bifurcation points
abstract type AbstractBifurcationPoint{U<:Real} end

# Sub-hierarchy for Limit Points
abstract type AbstractLimitPoint{U<:Real} <: AbstractBifurcationPoint{U} end

"""
    LimitPoint{U<:Real}

Generic limit point of an equilibrium branch. Characterized by a solution vector `x`,
the Jacobian matrix, the null vector of the extended system, and its index within the 
equilibrium branch.
"""
struct LimitPoint{U<:Real} <: AbstractLimitPoint{U}
    x::Vector{U}                  # Coordinates
    jacobian::Matrix{U}           # Jacobian matrix
    nullvec::Vector{U}            # Null vector of the extended system
    idx::Int                      # Index in the equilibrium branch
end

"""
    HamiltonianLimitPoint{U<:Real}

Limit point for a Hamiltonian system along an equilibrium branch. Includes energy `E` 
in addition to the regular fields.
"""
struct HamiltonianLimitPoint{U<:Real} <: AbstractLimitPoint{U}
    x::Vector{U}
    E::U                          # Energy at bifurcation
    jacobian::Matrix{U}
    nullvec::Vector{U}
    tol::U
    idx::Int
end

# Sub-hierarchy for Branch Points
abstract type AbstractBranchPoint{U<:Real} <: AbstractBifurcationPoint{U} end

"""
    BranchPoint{U<:Real}

Generic branch point where the equilibrium branch bifurcates. Includes a direction 
vector `dir`in the original equilibrium branch.
"""
struct BranchPoint{U<:Real} <: AbstractBranchPoint{U}
    x::Vector{U}
    jacobian::Matrix{U}
    nullvec::Vector{U}
    dir::Vector{U}                # Direction in the original equilibrium branch
    idx::Int
end

"""
    HamiltonianBranchPoint{U<:Real}

Branch point in a Hamiltonian context, includes the bifurcation energy `E`, along an 
equilibrium branch.
"""
struct HamiltonianBranchPoint{U<:Real} <: AbstractBranchPoint{U}
    x::Vector{U}
    E::U
    jacobian::Matrix{U}
    nullvec::Vector{U}
    dir::Vector{U}
    tol::U
    idx::Int
end

"""
    HopfPoint{U<:Real}

Hopf bifurcation point on an equilibrium branch. The null vector corresponds to a 
complex eigenvector whose associated eigenvalue has zero real part and nonzero imaginary 
part.
"""
struct HopfPoint{U<:Real} <: AbstractBifurcationPoint{U}
    x::Vector{U}
    jacobian::Matrix{U}
    nullvec::Vector{Complex{U}}  # Complex eigenvector associated with purely imaginary eigenvalue
    idx::Int
end


# -----------------------------
# Visual hierarchy (textual)
# -----------------------------

# AbstractBifurcationPoint{U}
# ├── AbstractLimitPoint{U}
# │   ├── LimitPoint{U}
# │   └── HamiltonianLimitPoint{U}
# ├── AbstractBranchPoint{U}
# │   ├── BranchPoint{U}
# │   └── HamiltonianBranchPoint{U}
# └── HopfPoint{U}