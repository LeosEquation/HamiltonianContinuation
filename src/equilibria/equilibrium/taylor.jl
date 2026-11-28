function update_taylor!(::EquilibriumProblem, workspace, problem)

    (; taylor, state) = workspace
    (; f!, params, t0, dim) = problem
    (; x_jet_1, dx_jet_1) = taylor
    (; x) = state

    for i in 1:dim
        x_jet_1[i][0][1] = x[i]
    end

    f!(dx_jet_1, x_jet_1, params, t0)

    return nothing

end