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
    state::ContinuationState,
    options::ContinuationOptions{U}
) where {U<:Real}

    (; dim) = problem
    (; maxsteps) = options

    branch = EquilibriumBranch(
        Array{U}(undef, dim, maxsteps),
        Array{Complex{U}}(undef, dim, maxsteps),
        Vector{Bool}(undef, maxsteps),
        Vector{U}(undef, maxsteps),
        LimitPoint{U}[],
        BranchPoint{U}[]
    )

    # Actualizar el primer paso con state

    return branch

end