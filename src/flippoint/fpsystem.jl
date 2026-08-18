function FPSystem!(
    FFP::Array{U,1},
    xTFP::Array{TaylorN{U},2},
    xFP::Array{U,1},
    x0FP::Array{U,1},
    dx0FP::Array{U,1},
    Δτ::U,
    nξ::T,
    nξFP::T,
    m::T,
    NFP::T,
) where {U<:Real,T<:Integer}

    comp3 = 0
    for comp1 in 1:m-1
        for comp2 in 1:nξFP
            comp3 += 1
            FFP[comp3] = xTFP[comp1, comp2][0][1] - xFP[nξFP+comp3]
        end
    end

    comp4 = 0
    for comp2 in 1:nξ
        comp3 += 1
        comp4 += 1
        FFP[comp3] = xFP[comp3] - xFP[comp4]
    end
    for comp2 in 1:nξ
        comp3 += 1
        comp4 += 1
        FFP[comp3] = xFP[comp3] + xFP[comp4]
    end

    FFP[NFP-1] = 0.0 # sum((xFP[i] - x0FP[i]) * dx0FP[i] for i in 1:nξ)
    FFP[NFP] = -1.0

    comp2 = 0

    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP-1] += Δτ * 0.5 * xFP[comp2] * dx0FP[comp2]
    end
    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP] += Δτ * 0.5 * xFP[comp2]^2
    end

    for comp3 in 2:m-1

        for comp1 in 1:nξ
            comp2 += 1
            FFP[NFP-1] += Δτ * xFP[comp2] * dx0FP[comp2]
        end
        for comp1 in 1:nξ
            comp2 += 1
            FFP[NFP] += Δτ * xFP[comp2]^2
        end

    end

    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP-1] += Δτ * 0.5 * xFP[comp2] * dx0FP[comp2]
    end
    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP] += Δτ * 0.5 * xFP[comp2]^2
    end

    return nothing
end

function FPJacobian!(
    JFP::AbstractArray{U,2},
    xTFP::Array{TaylorN{U},2},
    xFP::Array{U,1},
    dx0FP::Array{U,1},
    Δτ::U,
    nξ::T,
    nξFP::T,
    nλ::T,
    m::T,
    NξFP::T,
    NFP::T,
) where {U<:Real,T<:Integer}

    for point in 1:(m-1)

        ini = (point - 1) * nξFP

        for row in 1:nξFP

            for column in 1:nξFP
                JFP[ini+row, ini+column] = xTFP[point, row][1][column]
            end

            JFP[ini+row, ini+nξFP+row] = -one(U)

            for column in 1:nλ
                JFP[ini+row, NξFP+column] = xTFP[point, row][1][nξFP+column]
            end

        end

    end

    ### Set the coupling block

    ini = (m - 1) * nξFP

    for id in 1:nξ
        JFP[ini+id, ini+id] = one(U)
        JFP[ini+id, id] = -one(U)
    end

    for id in (nξ+1):nξFP
        JFP[ini+id, ini+id] = one(U)
        JFP[ini+id, id] = one(U)
    end

    ### Set the last two rows

    comp1 = 0

    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP-1, comp1] = 0.5 * Δτ * dx0FP[comp1]
    end
    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP, comp1] = Δτ * xFP[comp1]
    end

    for comp3 in 2:m-1

        for comp2 in 1:nξ
            comp1 += 1
            JFP[NFP-1, comp1] = Δτ * dx0FP[comp1]
        end
        for comp2 in 1:nξ
            comp1 += 1
            JFP[NFP, comp1] = 2.0 * Δτ * xFP[comp1]
        end

    end

    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP-1, comp1] = 0.5 * Δτ * dx0FP[comp1]
    end
    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP, comp1] = Δτ * xFP[comp1]
    end

    return nothing

end


function FPSystem!(
    FFP::Array{U,1},
    xTFP::Array{TaylorN{U},2},
    xFP::Array{U,1},
    x0FP::Array{U,1},
    dx0FP::Array{U,1},
    Φ::Array{U,1},
    Δτ::U,
    Δs::U,
    nξ::T,
    nξFP::T,
    NξFP::T,
    m::T,
    NFP::T,
) where {U<:Real,T<:Integer}

    comp3 = 0
    for comp1 in 1:m-1
        for comp2 in 1:nξFP
            comp3 += 1
            FFP[comp3] = xTFP[comp1, comp2][0][1] - xFP[nξFP+comp3]
        end
    end

    comp4 = 0
    for comp2 in 1:nξ
        comp3 += 1
        comp4 += 1
        FFP[comp3] = xFP[comp3] - xFP[comp4]
    end
    for comp2 in 1:nξ
        comp3 += 1
        comp4 += 1
        FFP[comp3] = xFP[comp3] + xFP[comp4]
    end

    FFP[NFP-2] = 0.0
    FFP[NFP-1] = -1.0

    comp2 = 0

    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP-2] += Δτ * 0.5 * xFP[comp2] * dx0FP[comp2]
    end
    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP-1] += Δτ * 0.5 * xFP[comp2]^2
    end

    for comp3 in 2:m-1

        for comp1 in 1:nξ
            comp2 += 1
            FFP[NFP-2] += Δτ * xFP[comp2] * dx0FP[comp2]
        end
        for comp1 in 1:nξ
            comp2 += 1
            FFP[NFP-1] += Δτ * xFP[comp2]^2
        end

    end

    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP-2] += Δτ * 0.5 * xFP[comp2] * dx0FP[comp2]
    end
    for comp1 in 1:nξ
        comp2 += 1
        FFP[NFP-1] += Δτ * 0.5 * xFP[comp2]^2
    end

    I = Δτ * (0.5 * sum((xFP[i] - x0FP[i]) * Φ[i] for i in 1:nξFP) +
              sum((xFP[i] - x0FP[i]) * Φ[i] for i in (nξFP+1):(NξFP-nξFP)) +
              0.5 * sum((xFP[i] - x0FP[i]) * Φ[i] for i in (NξFP-nξFP+1):(NξFP))) +
        (xFP[NFP-2] - x0FP[NFP-2]) * Φ[NFP-2] +
        (xFP[NFP-1] - x0FP[NFP-1]) * Φ[NFP-1] +
        (xFP[NFP] - x0FP[NFP]) * Φ[NFP]

    FFP[NFP] = I - Δs

    return nothing
end

function FPJacobian!(
    JFP::AbstractArray{U,2},
    xTFP::Array{TaylorN{U},2},
    xFP::Array{U,1},
    dx0FP::Array{U,1},
    Φ::Array{U,1},
    Δτ::U,
    nξ::T,
    nξFP::T,
    nλ::T,
    m::T,
    NξFP::T,
    NFP::T,
) where {U<:Real,T<:Integer}

    for point in 1:(m-1)
        ini = (point - 1) * nξFP
        for row in 1:nξFP
            for column in 1:nξFP
                JFP[ini+row, ini+column] = xTFP[point, row][1][column]
                JFP[ini+row, ini+nξFP+row] = -one(U)
            end
            for column in 1:nλ
                JFP[ini+row, NξFP+column] = xTFP[point, row][1][nξFP+column]
            end
        end
    end

    ### Set the coupling block

    ini = (m - 1) * nξFP

    for id in 1:nξ
        JFP[ini+id, ini+id] = one(U)
        JFP[ini+id, id] = -one(U)
    end

    for id in (nξ+1):nξFP
        JFP[ini+id, ini+id] = one(U)
        JFP[ini+id, id] = one(U)
    end

    ### Set the last two rows

    comp1 = 0

    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP-2, comp1] = 0.5 * Δτ * dx0FP[comp1]
    end
    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP-1, comp1] = Δτ * xFP[comp1]
    end

    for comp3 in 2:m-1

        for comp2 in 1:nξ
            comp1 += 1
            JFP[NFP-2, comp1] = Δτ * dx0FP[comp1]
        end
        for comp2 in 1:nξ
            comp1 += 1
            JFP[NFP-1, comp1] = 2.0 * Δτ * xFP[comp1]
        end

    end

    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP-2, comp1] = 0.5 * Δτ * dx0FP[comp1]
    end
    for comp2 in 1:nξ
        comp1 += 1
        JFP[NFP-1, comp1] = Δτ * xFP[comp1]
    end

    for i in 1:nξFP
        JFP[NFP, i] = 0.5 * Δτ * Φ[i]
    end

    for i in (nξFP+1):(NξFP-nξFP)
        JFP[NFP, i] = Δτ * Φ[i]
    end

    for i in (NξFP-nξFP+1):(NξFP)
        JFP[NFP, i] = 0.5 * Δτ * Φ[i]
    end

    JFP[NFP, NFP-2] = Φ[NFP-2]
    JFP[NFP, NFP-1] = Φ[NFP-1]
    JFP[NFP, NFP] = Φ[NFP]

    return nothing

end







function PDSystem!(FPD::Array{U,1},
    yT::Array{TaylorN{U},2},
    y::Array{U,1},
    y0::Array{U,1},
    dy0::Array{U,1},
    Φ0::Array{U,1},
    Δs::U,
    Δτ::U,
    n::T, m::T, ny::T, Ny::T) where {U<:Real,T<:Integer}

    comp3 = 0
    for comp1 in 1:m-1
        for comp2 in 1:ny-3
            comp3 += 1
            FPD[comp3] = yT[comp1, comp2][0][1] - y[ny-3+comp3]
        end
    end

    comp4 = 0
    for comp2 in 1:n-3
        comp3 += 1
        comp4 += 1
        FPD[comp3] = y[comp3] - y[comp4]
    end
    for comp2 in 1:n-3
        comp3 += 1
        comp4 += 1
        FPD[comp3] = y[comp3] + y[comp4]
    end

    FPD[Ny-2] = sum((y[i] - y0[i]) * dy0[i] for i in 1:n-3)
    FPD[Ny-1] = sum(y[i]^2 for i in n-2:ny-3) - 1.0
    FPD[Ny] = Δτ * sum((y[i] - y0[i]) * Φ0[i] for i in 1:ny-3) + sum((y[i] - y0[i]) * Φ0[i] for i in (ny-2):Ny) - Δs


    return nothing
end

function PDJacobian!(JPD::AbstractArray{U,2},
    yT::Array{TaylorN{U},2},
    y::Array{U,1},
    dy0::Array{U,1},
    Φ0::Array{U,1},
    Δτ::U,
    n::T, m::T, ny::T, Ny::T) where {U<:Real,T<:Integer}

    n_minus_2 = ny - 3
    oneU = one(U)

    # @show n

    @inbounds for comp1 in 1:m
        comp1_offset = (comp1 - 1) * n_minus_2
        for comp2 in 1:n_minus_2
            comp2_idx = comp1_offset + comp2

            if comp1 < m
                # Set the main block
                for comp3 in 1:n_minus_2
                    JPD[comp2_idx, comp1_offset+comp3] = yT[comp1, comp2][1][comp3]
                end

                JPD[comp2_idx, comp2_idx+n_minus_2] = -oneU

                JPD[comp2_idx, Ny-2] = yT[comp1, comp2][1][ny-2]
                JPD[comp2_idx, Ny-1] = yT[comp1, comp2][1][ny-1]
                JPD[comp2_idx, Ny] = yT[comp1, comp2][1][ny]
            else
                # @show comp2_idx
                # @show comp2
                if comp2 <= n - 3
                    JPD[comp2_idx, comp2] = -oneU
                else
                    JPD[comp2_idx, comp2] = oneU
                end
                JPD[comp2_idx, comp2_idx] = oneU
            end

        end
    end

    # Set the second-to-last row
    @inbounds for i in 1:n-3
        JPD[Ny-2, i] = dy0[i]
    end

    # Set the last row
    @inbounds for i in n-2:ny-3
        JPD[Ny-1, i] = 2 * y[i]
    end

    @inbounds for i in 1:ny-3
        JPD[Ny, i] = Δτ * Φ0[i]
    end

    @inbounds for i in (ny-2):Ny
        JPD[Ny, i] = Φ0[i]
    end


    return nothing
end