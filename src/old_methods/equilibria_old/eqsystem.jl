function EquilibriaJacobian!(J::Array{U, 2}, dx::Array{TaylorN{U}, 1}, 
                             Φ::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    for j in 1:n
        for i in 1:n-1
            @inbounds J[i,j] = dx[i][1][j]
        end
    end
    @inbounds J[n, :] .= Φ

    # GC.gc()

    return nothing

end

function EquilibriaSystem!(F::Array{U, 1}, dx::Array{TaylorN{U}, 1},
                           x1::Array{U, 1}, x0::Array{U, 1}, Φ::Array{U, 1},
                           Δs::U, n::T) where {U<:Real, T<:Integer}

    for i in 1:n-1
        @inbounds F[i] = dx[i][0][1]
    end
    @inbounds F[n] = sum((x1[i] - x0[i])*Φ[i] for i in 1:n) - Δs

    # GC.gc()

    return nothing

end

function HEquilibriaJacobian!(J::Array{U, 2}, H::TaylorN{U}, dx::Array{TaylorN{U}, 1}, 
                             Φ::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    for j in 1:n
        for i in 1:n-2
            @inbounds J[i,j] = dx[i][1][j]
        end
        J[n-1, j] = H[1][j]
    end
    
    @inbounds J[n, :] .= Φ

    # GC.gc()

    return nothing

end

function HEquilibriaSystem!(F::Array{U, 1}, Heval::U, E::U, dx::Array{TaylorN{U}, 1},
                           x1::Array{U, 1}, x0::Array{U, 1}, Φ::Array{U, 1},
                           Δs::U, n::T) where {U<:Real, T<:Integer}

    for i in 1:n-2
        @inbounds F[i] = dx[i][0][1]
    end
    F[n-1] = Heval - E
    @inbounds F[n] = sum((x1[i] - x0[i])*Φ[i] for i in 1:n) - Δs

    # GC.gc()

    return nothing

end