abstract type EquilibriumBifurcations{U<:Real} end

struct LimitPoint{U<:Real} <: EquilibriumBifurcations

    x::Vector{U}
    eigenvalues::Vector{Complex{U}}
    tangent::Vector{U}

end

struct BranchPoint{U<:Real} <: EquilibriumBifurcations

    x::Vector{U}
    eigenvalues::Vector{Complex{U}}
    tangent::Vector{U}
    emerging_tangent::Vector{U}

end

abstract type AbstractBifurcationTest{U<:Real} end

struct LimitPointTest{U<:Real} <: AbstractBifurcationTest
    change::Vector{U}
end

struct BranchPointTest{U<:Real} <: AbstractBifurcationTest
    change::Vector{U}
end

function build_test()
    return [
        LimitPointTest{U}(zeros(U, 2)),
        BranchPointTest{U}(zeros(U, 2))
    ]
end