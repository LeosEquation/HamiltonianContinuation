function EquilibriaFinding!(f!, x_ini::Array{U,1}, params, tol::U, ite::T) where {U<:Real,T<:Integer}

    n = length(x_ini)

    δx = TaylorN.(1:n, order=1)

    Δx = zeros(U, n)

    xaux = copy(δx)
    dx = zero(δx)

    Jeval = zeros(U, n, n)
    dxeval = zeros(U, n)

    for comp1 in 1:n
        xaux[comp1][0][1] = x_ini[comp1]
    end

    f!(dx, xaux, params, 0.0)

    for comp1 in 1:n
        for comp2 in 1:n
            Jeval[comp1, comp2] = dx[comp1][1][comp2]
        end
        dxeval[comp1] = dx[comp1][0][1]
    end

    j = 1

    while j <= ite && norm(dxeval) > tol

        Δx = Jeval \ dxeval

        for comp1 in 1:n
            x_ini[comp1] -= Δx[comp1]
            xaux[comp1][0][1] = x_ini[comp1]
        end

        f!(dx, xaux, params, 0.0)

        for comp1 in 1:n
            for comp2 in 1:n
                Jeval[comp1, comp2] = dx[comp1][1][comp2]
            end
            dxeval[comp1] = dx[comp1][0][1]
        end

        j += 1

    end

    return nothing

end