function build_variables(
    ::Type{EquilibriumContinuation},
    problem::EquilibriumProblem
)
    return copy(problem.x_ini)
end