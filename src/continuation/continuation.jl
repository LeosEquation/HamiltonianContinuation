"""
    Continuation(::Type{C}, problem, options)

Performs a continuation process using the specified continuation method.

Initializes the continuation workspace and iteratively computes new
continuation points using a predictor-corrector scheme with Newton's
method.

After each accepted step, updates the tangent vector, stores the new
continuation point, detects bifurcations, and updates the previous
continuation state. The continuation step size is adapted according
to the number of Newton iterations required for convergence.

The process terminates when the maximum number of continuation steps
is reached or a stopping condition is met.

# Arguments
- `C::Type{<:AbstractContinuation}`: Continuation method to be used.
- `problem::AbstractProblem`: Problem defining the system to be continued.
- `options::ContinuationOptions`: Options controlling the continuation
  algorithm.

# Returns
The computed continuation branch, containing the accepted continuation
points up to the last completed step.
"""
function Continuation(
    ::Type{C},
    problem::P,
    space::JetSpace,
    options::ContinuationOptions{U}
) where {C<:AbstractContinuation,P<:AbstractProblem,U<:Real}

    # Initialization

    workspace, Δs_current, maxsteps =
        build_continuation_workspace(C, space, problem, options)

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