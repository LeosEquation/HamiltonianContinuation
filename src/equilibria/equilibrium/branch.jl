struct EquilibriumBranch{U<:Real}
    x::Array{U,2}
    values::Array{Complex{U},2}
    stability::Vector{Bool}
    residual::Array{U,1}
    limit_point::Vector{LimitPoint{U}}
    branch_point::Vector{BranchPoint{U}}
end

function build_branch(
    ::Type{EquilibriumContinuation},
    problem::EquilibriumProblem{U},
    state::ContinuationState{U},
    phase_space::EquilibriumPhaseSpace{U},
    options::ContinuationOptions{U}
) where {U<:Real}

    (; dim) = problem
    (; maxsteps) = options

    branch = EquilibriumBranch(
        Matrix{U}(undef, dim, maxsteps),
        Matrix{Complex{U}}(undef, dim, maxsteps),
        Vector{Bool}(undef, maxsteps),
        Vector{U}(undef, maxsteps),
        LimitPoint{U}[],
        BranchPoint{U}[]
    )

    # Store the initial continuation point
    branch.x[:, 1] .= state.x
    branch.values[:, 1] .= eigvals(phase_space.J)
    branch.residual[1] = norm(state.F, Inf)

    (; atol, rtol) = options

    λ = eigvals(phase_space.J)
    tol = atol + rtol * norm(state.x, Inf)

    branch.stability[1] = !any(real.(λ) .> tol)

    return branch
end