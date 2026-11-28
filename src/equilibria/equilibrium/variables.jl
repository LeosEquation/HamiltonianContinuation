function build_variables(
    ::Type{EquilibriumContinuation},
    problem::EquilibriumProblem{U}
) where {U<:Real}

    return copy(problem.x_ini)

end