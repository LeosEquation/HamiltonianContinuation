struct ContinuationWorkspace{B,S,F}
    branch::B
    state::S
    bifurcations::F
end

function build_continuation_workspace(
    ::Type{C},
    problem,
    options
) where {C<:AbstractContinuation}

    branch = build_branch(C, problem, options)

    state = build_continuation_state(C, problem)

    bifurcations = build_bifurcations(C, problem)

    Δs_current = options.Δs

    maxsteps = options.maxsteps

    return ContinuationWorkspace(
        branch,
        state,
        bifurcations
    ), Δs_current, maxsteps

end