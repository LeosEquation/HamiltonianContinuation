function build_jacobian(
    ::EquilibriumContinuation,
    taylor,
    problem::EquilibriumProblem{U}
) where {U<:Real}

    (; dim_variables, dim) = problem
    (; dx_jet_1) = taylor

    n = dim_variables

    # Allocate the extended Jacobian
    J = zeros(U, dim, dim)

    # Compute the physical Jacobian and parameter derivatives
    @inbounds for j in 1:dim
        for i in 1:n
            J[i, j] = dx_jet_1[i][1][j]
        end
    end

    # Compute the null space of the initial Jacobian
    T0 = nullspace(J)

    size(T0, 2) == 1 ||
        throw(ArgumentError(
            "The initial Jacobian must have a one-dimensional null space."
        ))

    # Extract the null vector
    T0 = vec(T0)

    # Complete the extended Jacobian
    @inbounds J[end, :] .= T0

    return J, T0
end