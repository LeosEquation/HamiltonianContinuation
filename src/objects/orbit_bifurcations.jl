# Complete hierarchy for bifurcation points in general and Hamiltonian systems

# Sub-hierarchy for Period Doubling
abstract type AbstractFlipPoint{U<:Real} <: AbstractBifurcationPoint{U} end

"""
    FlipPoint{U<:Real}

Generic flip (period-doubling) point in a branch of periodic orbits.
Characterized by:
  • The time partition `tⱼ = t[j]`, normalized over the interval [0, 1]
  • The discrete solution (including the bifurcation parameter and the period of each orbit)
    `xᵢ(tⱼ) = x[j, i]`
  • The corresponding discrete eigenfunction `vᵢ(tⱼ) = v[j, i]`
  • The index of the preceding point in the branch `idx`
"""
struct FlipPoint{U<:Real} <: AbstractFlipPoint{U}
    t::LinRange{U}
    x::Array{U,2}
    v::Array{U,2}
    idx::Int
    err::U
end

"""
    HamiltonianFlipPoint{U<:Real}

Generic flip (period-doubling) point in a branch of Hamiltonian periodic orbits.
Characterized by:
    - The time partition `tⱼ = t[j]`, normalized over the interval [0, 1]
    - The discrete solution (including the bifurcation parameter and the period of each orbit)
        `xᵢ(tⱼ) = x[j, i]`
    - The constant energy `E` associated with the orbit
    - The corresponding discrete eigenfunction `vᵢ(tⱼ) = v[j, i]`
    - The index of the preceding point in the branch `idx`
"""
struct HamiltonianFlipPoint{U<:Real} <: AbstractFlipPoint{U}
    t::LinRange{U}
    x::Array{U,2}
    E::U
    v::Array{U,2}
    μ::Array{Complex{U},1}
    idx::Int
    tol::U
end


# Sub-hierarchy for Period Doubling
abstract type AbstractBranchCyclePoint{U<:Real} <: AbstractBifurcationPoint{U} end

"""
    LimitCyclePoint{U<:Real}

Generic limit cycle point in a branch of periodic orbits.
Characterized by:
  • The time partition `tⱼ = t[j]`, normalized over the interval [0, 1]
  • The discrete solution (including the bifurcation parameter and the period of each orbit)
    `xᵢ(tⱼ) = x[j, i]`
  • The corresponding discrete eigenfunction `vᵢ(tⱼ) = v[j, i]`
  • The index of the preceding point in the branch `idx`
"""
struct BranchCyclePoint{U<:Real} <: AbstractBranchCyclePoint{U}
    t::LinRange{U}
    x::Array{U,2}
    v::Array{U,2}
    idx::Int
    tol::U
end

"""
    HamiltonianBranchCyclePoint{U<:Real}

Generic flip (period-doubling) point in a branch of Hamiltonian periodic orbits.
Characterized by:
    - The time partition `tⱼ = t[j]`, normalized over the interval [0, 1]
    - The discrete solution (including the bifurcation parameter and the period of each orbit)
        `xᵢ(tⱼ) = x[j, i]`
    - The constant energy `E` associated with the orbit
    - The corresponding discrete eigenfunction `vᵢ(tⱼ) = v[j, i]`
    - The index of the preceding point in the branch `idx`
"""
struct HamiltonianBranchCyclePoint{U<:Real} <: AbstractBranchCyclePoint{U}
    t::LinRange{U}
    x::Array{U,2}
    E::U
    v::Array{U,2}
    idx::Int
    tol::U
end

