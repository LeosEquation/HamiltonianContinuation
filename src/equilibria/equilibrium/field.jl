function build_field(
    ::EquilibriumContinuation,
    taylor,
    problem::EquilibriumProblem{U}
) where {U<:Real}

    (; dim_variables, dim) = problem
    (; dx_jet_1) = taylor

    F = zeros(U, dim)

    # Evaluate the physical vector field
    @inbounds for i in 1:dim_variables
        F[i] = dx_jet_1[i][0][1]
    end

    # Initial continuation condition
    F[dim] = zero(U)

    return F
end

# Quiero que sea el análogo a esta función:
# function EquilibriaSystem!(F::Array{U, 1}, dx::Array{TaylorN{U}, 1},
#                            x1::Array{U, 1}, x0::Array{U, 1}, Φ::Array{U, 1},
#                            Δs::U, n::T) where {U<:Real, T<:Integer}

#     for i in 1:n-1
#         @inbounds F[i] = dx[i][0][1]
#     end
#     @inbounds F[n] = sum((x1[i] - x0[i])*Φ[i] for i in 1:n) - Δs

#     # GC.gc()

#     return nothing

# end