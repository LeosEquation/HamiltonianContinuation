
function LPJacobian!(J::Array{U, 2}, Jeval::Array{U, 2}, dx::Array{TaylorN{U}, 1},
                     q::Array{U, 1}, Φ::Array{U, 1}, 
                     ordtup::Array{NTuple{N, T}, 2}, n::T) where {U<:Real, T<:Integer, N}

    i1 = n - 2
    i2 = 2*n - 3

    for i in 1:n-2

        for j in 1:n
            J[i, j] = Jeval[i, j]
            # J[i1+i, j] = sum(differentiate(ntuple(l -> count(==(l), (k, j)), TaylorSeries.get_numvars()), dx[i]) * q[n+k] for k in 1:n-2)
            J[i1+i, j] = sum(differentiate(ordtup[j, k], dx[i]) * q[n+k] for k in 1:n-2)
            if j < (n - 1) 
                J[i1+i, n+j] = Jeval[i, j]
            end
        end

        J[i2, n+i] = 2.0 * q[n+i]
    end

    for i in 1:n
        J[2*n-2, :] .= Φ
    end

    # GC.gc()

    return nothing

end

function LPSystem!(F::Array{U, 1}, Jeval::Array{U, 2}, dxeval::Array{U, 1},
                   q::Array{U, 1}, q0::Array{U, 1}, Φ::Array{U, 1}, Δs::U, n::T) where {U<:Real, T<:Integer}

    i1 = n-2
    for i in 1:n-2
        F[i] = dxeval[i]
        F[i1+i] = sum(Jeval[i, k] * q[n+k] for k in 1:n-2)
    end
    F[2*n-3] = sum(q[n+i]^2 for i in 1:n-2) - 1.0
    F[2*n-2] = sum((q[i] - q0[i]) * Φ[i] for i in 1:(2*n-2)) - Δs

    # GC.gc()

    return nothing

end

function LPJacobian!(J::Array{U, 2}, Jeval::Array{U, 2}, dx::Array{TaylorN{U}, 1},
                     q::Array{U, 1}, ordtup::Array{NTuple{N, T}, 2}, n::T) where {U<:Real, T<:Integer, N}

    i1 = n-1
    i2 = 2*n - 1

    for i in 1:n-1
        for j in 1:n
            J[i, j] = Jeval[i, j]
            # J[i1+i, j] = sum(differentiate(ntuple(l -> count(==(l), (k, j)), TaylorSeries.get_numvars()), dx[i]) * q[n+k] for k in 1:n-1)
            J[i1+i, j] = sum(differentiate(ordtup[j, k], dx[i]) * q[n+k] for k in 1:n-1)
            if j < n
                J[i1+i, n+j] = Jeval[i, j]
            end
        end
        J[i2, n+i] = 2.0 * q[n+i]
    end

    # GC.gc()

    return nothing

end

function LPSystem!(F::Array{U, 1}, Jeval::Array{U, 2}, dxeval::Array{U, 1},
                   q::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    i1 = n-1
    for i in 1:n-1
        F[i] = dxeval[i]
        F[i1+i] = sum(Jeval[i, k] * q[n+k] for k in 1:n-1)
    end
    F[2*n-1] = sum(q[n+i]^2 for i in 1:n-1) - 1.0

    # GC.gc()

    return nothing

end