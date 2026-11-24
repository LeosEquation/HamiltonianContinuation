function OrbitSystem!(  
                        F::Array{U, 1},
                        xT::Array{TaylorN{U}, 2},
                        x1::Array{U, 1},
                        x0::Array{U, 1},
                        dx0::Array{U, 1},
                        Φ::Array{U, 1},
                        Δs::U,
                        Δτ::U,
                        nξ::T,
                        m::T,
                        Nξ::T,
                        N::T
                        ) where {U<:Real, T<:Integer}

    comp3 = 0        
    for comp1 in 1:m-1
        for comp2 in 1:nξ
            comp3 += 1
            F[comp3] = xT[comp1, comp2][0][1] - x1[nξ + comp3]
        end
    end

    for comp2 in 1:nξ
        comp3 += 1
        F[comp3] = x1[comp3] - x1[comp2]
    end

    F[N-1] = Δτ * ( 0.5 * sum( x1[i] * dx0[i] for i in 1:nξ ) +
                    sum( x1[i] * dx0[i] for i in (nξ+1):(Nξ-nξ) ) +
                    0.5 * sum( x1[i] * dx0[i] for i in (Nξ-nξ+1):(Nξ) ) )

    # F[N-1] = sum((x1[i] - x0[i]) * dx0[i] for i in 1:nξ)

    I = Δτ * ( 0.5 * sum( (x1[i] - x0[i]) * Φ[i] for i in 1:nξ ) +
               sum( (x1[i] - x0[i]) * Φ[i] for i in (nξ+1):(Nξ-nξ) ) +
               0.5 * sum( (x1[i] - x0[i]) * Φ[i] for i in (Nξ-nξ+1):(Nξ) ) )
        + (x1[N-1] - x0[N-1]) * Φ[N-1] + (x1[N] - x0[N]) * Φ[N]

    F[N]   = I - Δs

    return nothing

end

function OrbitJacobian!(
                        J::AbstractArray{U,2}, 
                        xT::Array{TaylorN{U},2}, 
                        dx0::Array{U,1}, 
                        Φ::Array{U,1}, 
                        Δτ::U, 
                        nξ::T, 
                        nλ::T,
                        m::T,
                        Nξ::T, 
                        N::T
                        ) where {U<:Real, T<:Integer}

    ### Set the main blocks

    for point in 1:(m-1)
        ini = (point-1)*nξ
        for row in 1:nξ
            for column in 1:nξ
                J[ini + row, ini + column] = xT[point, row][1][column]
                J[ini + row, ini + nξ + row] = - one(U) 
            end
            for column in 1:nλ
                J[ini + row, Nξ + column] = xT[point, row][1][nξ + column]
            end
        end
    end

    ### Set the coupling block

    ini = (m-1)*nξ
    for id in 1:nξ
        J[ini + id, ini + id] = one(U)
        J[ini + id, id] = - one(U)
    end

    ### Set the penultimate row

    @inbounds for i in 1:nξ
        J[N-1, i] = 0.5 * Δτ * dx0[i]
    end

    @inbounds for i in (nξ+1):(Nξ-nξ)
        J[N-1, i] = Δτ * dx0[i]
    end

    @inbounds for i in (Nξ-nξ+1):(Nξ)
        J[N-1, i] = 0.5 * Δτ * dx0[i]
    end

    ### Set the last row

    @inbounds for i in 1:nξ
        J[N, i] = 0.5 * Δτ * Φ[i]
    end

    @inbounds for i in (nξ+1):(Nξ-nξ)
        J[N, i] = Δτ * Φ[i]
    end

    @inbounds for i in (Nξ-nξ+1):(Nξ)
        J[N, i] = 0.5 * Δτ * Φ[i]
    end

    J[N, N-1] = Φ[N-1]
    J[N, N] = Φ[N]
    
    return nothing

end