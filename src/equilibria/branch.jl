struct EquilibriumBranch{U<:Real}
    coordinates::Array{U,2}
    eigenvalues::Array{Complex{U},2}
    bifurcations::Array{EquilibriumBifurcations{U},1}
    stability::Vector{Bool}
    residual::Array{U,1}
end