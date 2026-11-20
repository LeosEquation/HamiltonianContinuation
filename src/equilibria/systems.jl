build_eq_system(dimension) = function System!(continuation_point, Δs)

    (; F, dx, x, x0, tangent0) = continuation_point

    for i in 1:(dimension-1)
        F[i] = dx[i][0][1]
    end

    F[end] = sum((x[i][0][1] - x0[i]) * tangent0[i]
                 for i in eachindex(x)) - Δs

    return nothing

end

build_eq_jacobian(dimension) = function Jacobian!(continuation_point)

    (; J, dx, tangent0) = continuation_point

    for j in 1:dimension
        for i in 1:(dimension-1)
            J[i, j] = dx[i][1][j]
        end
    end

    J[dimension, :] .= tangent0

    return nothing

end

build_lp_system(dimension) = function LP_System!(bifurcation_point)

    (; F, dx, q) = bifurcation_point

    ncoord = dimension - 1

    for i in 1:ncoord

        F[i] = dx[i][0][1]

        F[ncoord+i] = sum(
            dx[i][1][k] * q[dimension+k]
            for k in 1:ncoord
        )

    end

    F[2*dimension-1] = sum(
        q[dimension+i]^2
        for i in 1:ncoord
    ) - 1.0

    return nothing

end

build_lp_jacobian(dimension, ordtup) = function LP_Jacobian!(bifurcation_point)

    (; J, dx, q) = bifurcation_point

    nxcoord = dimension - 1
    nq = 2*dimension - 1

    for i in 1:nxcoord
        for j in 1:dimension

            J[i, j] = dx[i][1][j]

            J[nxcoord+i, j] = sum(
                differentiate(ordtup[j, k], dx[i]) * q[dimension+k]
                for k in 1:nxcoord
            )

            if j < dimension
                J[nxcoord+i, dimension+j] = dx[i][1][j]
            end

        end

        J[nq, dimension+i] = 2.0 * q[dimension+i]
    end

    return nothing

end

build_bp_system(dimension) = function BP_System!(bifurcation_point)

    (; F, dx, q) = bifurcation_point

    ncoord = dimension-1

    for i in 1:ncoord

        F[i] = dx[i][0][1] + q[2*dimension] * q[dimension+i]

        F[dimension-1+i] = sum(
            q[dimension+k] * dx[k][1][i]
            for k in 1:ncoord
        )

    end

    F[2*dimension-1] = sum(
        q[dimension+k] * dx[k][1][dimension]
        for k in 1:ncoord
    )

    F[2*dimension] = sum(
        q[dimension+k]^2
        for k in 1:ncoord
    ) - 1.0

    return nothing

end

build_bp_jacobian(dimension, ordtup) = function BP_Jacobian!(bifurcation_point)

    (; J, dx, q) = bifurcation_point

    nxcoord = dimension-1
    nq = 2*dimension
    nqcoord = nq - 1

    for j in 1:dimension

        for i in 1:nxcoord

            J[i, j] = dx[i][1][j]

            J[nxcoord+i, j] = sum(
                differentiate(ordtup[i, j], dx[k]) * q[dimension+k]
                for k in 1:nxcoord
            )

            if j < dimension
                J[nxcoord+i, dimension+j] = dx[j][1][i]
            end

        end

        J[nqcoord, j] = sum(
            differentiate(ordtup[dimension, j], dx[k]) * q[dimension+k]
            for k in 1:nxcoord
        )

        if j < dimension
            J[nqcoord, dimension+j] = dx[j][1][dimension]
        end

    end

    for j in 1:nxcoord

        J[j, dimension+j] = q[nq]
        J[j, nq] = q[dimension+j]
        J[nq, dimension+j] = 2.0 * q[dimension+j]

    end

    return nothing

end