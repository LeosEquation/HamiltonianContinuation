struct EquilibriumContinuationPoint{U<:Real}
    dx::Vector{TaylorN{U}}
    x::Vector{TaylorN{U}}
    params::Any
    t::U
    tangent::Vector{U}
    F::Vector{U}
    J::Matrix{U}
    x0::Vector{U}
    dx0::Vector{U}
    tangent0::Vector{U}
end

function build_continuation_point(problem::EquilibriumProblem, Δs)

    (; dimension, space, init, f!, params, t, System!, Jacobian!) = problem

    U = eltype(init)

    dx = variables(space[1])
    x = variables(space[1])

    tangent = zeros(U, dimension)
    tangent0 = zeros(U, dimension)

    x0 = copy(init)
    dx0 = zeros(U, dimension)

    for i in eachindex(x)
        x[i][0][1] = init[i]
    end

    f!(dx, x, params, t)
    f!(dx0, x0, params, t)

    F = zeros(U, dimension)
    J = zeros(U, dimension, dimension)

    continuation_point = EquilibriumContinuationPoint(
        dx, x, params, t, tangent, F, J, x0, dx0, tangent0
    )

    System!(continuation_point, Δs)
    Jacobian!(continuation_point)

    NS = nullspace(J)

    if size(NS, 2) != 1
        throw(ArgumentError(
            "The initial equilibrium is not a regular point."
        ))
    end

    tangent0 .= NS[:, 1]
    tangent .= NS[:, 1]

    Jacobian!(continuation_point)

    return continuation_point
end