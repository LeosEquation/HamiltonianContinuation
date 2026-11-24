function build_field(C, problem)

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