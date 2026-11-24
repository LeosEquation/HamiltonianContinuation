function Continuation(
    ::Type{C},
    problem::AbstractProblem,
    options::ContinuationOptions{U}
) where {C<:AbstractContinuation,U<:Real}

    # Initialization

    workspace, Δs_current, maxsteps =
        build_continuation_workspace(C, problem, options)

    step = 2
    terminated = false

    # Continuation loop

    while step <= maxsteps && !terminated

        # Predictor-corrector iteration

        accepted, terminated, ite, Δs_current = continuation_step!(
            C, workspace, problem, options, Δs_current
        )

        if terminated || !accepted
            break
        end

        # Update tangent

        update_tangent!(workspace)

        # Store accepted point

        update_branch!(C, workspace, step)

        # Detect bifurcations

        find_bifurcations!(
            C, workspace, problem, options
        )

        # Update continuation state

        update_old_state!(workspace)

        step += 1

        # Adapt step size

        Δs_current = adapt_stepsize(Δs_current, ite, options)

    end

    return build_branch(C, workspace, step - 1)

end