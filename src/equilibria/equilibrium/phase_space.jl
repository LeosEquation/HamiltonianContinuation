struct EquilibriumPhaseSpace{U<:Real}
    x::Vector{U}
    p::Vector{U}
    dx::Vector{U}
    J::Matrix{U}
    dx_p::Matrix{U}
end

function build_phase_space(
    ::EquilibriumContinuation,
    taylor,
    problem::EquilibriumProblem{U}
) where {U<:Real}

    (; dim_variables) = problem
    (; x_jet_1, dx_jet_1) = taylor

    n = dim_variables

    # Allocate physical variables and parameters
    x = zeros(U, n)
    p = zeros(U, 1)

    # Allocate the vector field and its derivatives
    dx = similar(x)
    J = Matrix{U}(undef, n, n)
    dx_p = Matrix{U}(undef, n, 1)

    # Extract Taylor coefficients
    for i in 1:n
        x[i] = x_jet_1[i][0][1]
        dx[i] = dx_jet_1[i][0][1]
        dx_p[i, 1] = dx_jet_1[i][1][i]

        for j in 1:n
            J[i, j] = dx_jet_1[i][1][j]
        end
    end

    p[1] = x_jet_1[end][0][1]

    return EquilibriumPhaseSpace(x, p, dx, J, dx_p)
end

function update_phase_space!(
    ::EquilibriumContinuation,
    phase_space::EquilibriumPhaseSpace,
    state::ContinuationState
)
    (; x, p, dx, J, dx_p) = phase_space
    (; F) = state

    n = length(dx)

    copyto!(x, view(state.x, 1:n))
    p[1] = state.x[end]

    copyto!(dx, view(F, 1:n))
    copyto!(J, view(state.J, 1:n, 1:n))
    copyto!(dx_p, view(state.J, 1:n, n + 1))

    return nothing
end