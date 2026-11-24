function compute_correction!(
    workspace,
    Δs_current,
    problem,
    options
)

    update_field!(workspace, Δs_current, problem)
    update_jacobian!(workspace, problem)

    fail = not_regular(workspace, options)

    if !fail
        workspace.state.Δx .=
            workspace.state.J \ workspace.state.F
    end

    return fail

end

function continuation_step!(
    ::Type{C},
    workspace,
    problem,
    options,
    Δs_current
) where {C<:AbstractContinuation}

    accepted = false
    terminated = false
    retry = 0
    ite = 1

    while !accepted && retry <= options.max_retries

        # Prediction

        predict!(workspace, Δs_current)

        if stop(workspace, problem)
            terminated = true
            break
        end

        # Initial correction

        fail = compute_correction!(
            workspace, Δs_current, problem, options
        )

        ite = 1

        # Newton correction

        while !fail && ite < options.maxite &&
            !convergence(workspace, options)

            update_newton!(workspace)

            if stop(workspace, problem)
                terminated = true
                fail = true
                break
            end

            fail = compute_correction!(
                workspace, Δs_current, problem, options
            )

            if !fail
                ite += 1
            end

        end

        # Check convergence

        if !terminated && !fail &&
            convergence(workspace, options)

            accepted = true

        elseif !terminated

            Δs_current *= options.shrink_factor
            retry += 1

            if Δs_current < options.Δsmin
                terminated = true
            end

        end

    end

    return accepted, terminated, ite, Δs_current

end