
"""
    compute_correction!(::Type{C}, workspace, problem, Δs_current)

Computes the Newton correction for the current continuation point.

Updates the vector field and its Jacobian according to the continuation
method `C`, the problem definition, and the current continuation step size.
Then, computes an LU factorization of the Jacobian and solves the linear
system

    J * Δx = F

The factorization is performed on an auxiliary copy of the Jacobian,
preserving the original matrix. The resulting Newton correction is stored
in `Δx`.

# Arguments
- `C::Type{<:AbstractContinuation}`: Continuation method to be used.
- `workspace::ContinuationWorkspace`: Workspace containing the numerical
  state and auxiliary data.
- `problem::AbstractProblem`: Problem defining the system to be continued.
- `Δs_current::Real`: Current continuation step size.

# Returns
Returns `true` if the LU factorization is unsuccessful, indicating that
the Newton correction could not be computed, and `false` otherwise.
"""
function compute_correction!(
    ::Type{C},
    workspace::ContinuationWorkspace,
    problem::AbstractProblem,
    Δs_current::U
) where {C<:AbstractContinuation,U<:Real}

    (; state) = workspace
    (; Δx, J, F, J_factor) = state

    update_taylor!(C, workspace, problem)
    update_field!(C, workspace, Δs_current, problem)
    update_jacobian!(C, workspace, problem)

    copyto!(J_factor, J)
    luJ = lu!(J_factor; check=false)

    if !issuccess(luJ)
        return true
    end

    ldiv!(Δx, luJ, F)

    return false

end


"""
    predict!(workspace::ContinuationWorkspace, Δs_current)

Computes the predictor for the next continuation point using a
first-order extrapolation along the current tangent direction.

Updates the state vector according to

    x = x0 + Δs_current * T0

where `x0` is the previous continuation point and `T0` is its
associated tangent vector.

Returns `nothing`.
"""
function predict!(
    workspace::ContinuationWorkspace,
    Δs_current::U
) where {U<:Real}

    (; state) = workspace
    (; x, x0, T0) = state

    for i in eachindex(x)
        x[i] = x0[i] + Δs_current * T0[i]
    end

    return nothing

end


"""
    convergence(workspace::ContinuationWorkspace, options::ContinuationOptions)

Checks whether the current Newton iteration satisfies the convergence
criterion.

Convergence is determined using the infinity norm of the residual
vector `F` and a tolerance defined by

    tol = atol + rtol * norm(x, Inf)

where `atol` and `rtol` are the absolute and relative tolerances,
respectively, and `x` is the current state vector.

Returns `true` if `norm(F, Inf) <= tol`, and `false` otherwise.
"""
function convergence(
    workspace::ContinuationWorkspace,
    options::ContinuationOptions{U}
) where {U<:Real}
    (; state) = workspace
    (; F, x) = state
    (atol, rtol) = options

    tol = atol + rtol * norm(x, Inf)

    return norm(F, Inf) <= tol

end


"""
    continuation_step!(::Type{C}, workspace, problem, options, Δs_current)

Performs a single continuation step using predictor-corrector continuation
with Newton's method.

First, computes a predictor for the next continuation point. Then,
iteratively corrects the prediction until the convergence criterion is
satisfied or the correction procedure fails.

If the correction fails or the maximum number of Newton iterations is
reached without convergence, the continuation step is retried with a
reduced step size. The reduction factor and retry limits are determined
by `options`.

The process terminates if the continuation stopping condition is met,
the minimum step size is reached, or the retry limit is exceeded.

# Arguments
- `C::Type{<:AbstractContinuation}`: Continuation method to be used.
- `workspace::ContinuationWorkspace`: Workspace containing the branch,
  numerical state, and bifurcation data.
- `problem::AbstractProblem`: Problem defining the system to be continued.
- `options::ContinuationOptions`: Options controlling the continuation
  algorithm and Newton iterations.
- `Δs_current::Real`: Current continuation step size.

# Returns
A tuple containing:
- `accepted::Bool`: Indicates whether the continuation step converged
  successfully.
- `terminated::Bool`: Indicates whether the continuation process should
  terminate.
- `ite::Int`: Number of successfully computed Newton corrections in
  the final attempt.
- `Δs_current`: Continuation step size resulting from the accepted
  attempt or the final retry.

The workspace is updated during the procedure.
"""
function continuation_step!(
    ::Type{C},
    workspace::ContinuationWorkspace,
    problem::P,
    options::ContinuationOptions{U},
    Δs_current::U
) where {C<:AbstractContinuation,P<:AbstractProblem,U<:Real}

    (; state) = workspace
    (; maxite, max_retries, shrink_factor, Δsmin) = options
    (; x, Δx) = state

    accepted = false
    terminated = false
    retry = 0
    ite = 1

    while !accepted && retry <= max_retries

        # Prediction

        predict!(workspace, Δs_current)

        if stop(C, workspace, problem)
            terminated = true
            break
        end

        # Initial correction

        fail = compute_correction!(
            C, workspace, problem, options, Δs_current
        )

        ite = fail ? 0 : 1

        # Newton correction

        while !fail && ite < maxite &&
            !(convergence(workspace, options))

            x .-= Δx

            if stop(C, workspace, problem)
                terminated = true
                fail = true
                break
            end

            fail = compute_correction!(
                C, workspace, problem, options, Δs_current
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

            Δs_current *= shrink_factor
            retry += 1

            if Δs_current < Δsmin
                terminated = true
            end

        end

    end

    return accepted, terminated, ite, Δs_current

end

"""
    update_tangent!(workspace::ContinuationWorkspace)

Updates the tangent vector by solving the linear system

    J * T = v

using an LU factorization of the Jacobian.

The factorization is performed on an auxiliary copy of `J`, preserving
the original Jacobian. The resulting tangent vector is stored in `T`.

Returns `true` if the LU factorization is unsuccessful, and `false`
otherwise.
"""
function update_tangent!(workspace::ContinuationWorkspace)

    (; state) = workspace
    (; T, J, v, J_factor) = state

    copyto!(J_factor, J)
    luJ = lu!(J_factor; check=false)

    if !issuccess(luJ)
        return true
    end

    ldiv!(T, luJ, v)

    return false

end

"""
    update_old_state!(workspace::ContinuationWorkspace)

Updates the previous continuation state by copying the current
continuation point and tangent vector into their corresponding
previous-state arrays.

Specifically, updates `x0` with `x` and `T0` with `T`.

Returns `nothing`.
"""
function update_old_state!(workspace::ContinuationWorkspace)

    (; state) = workspace
    (; x, x0, T, T0) = state

    copyto!(x0, x)
    copyto!(T0, T)

    return nothing

end


"""
    adapt_stepsize(Δs_current, ite, options)

Adapts the continuation step size according to the number of Newton
iterations required to converge the current continuation step.

If the number of iterations is below the target, the step size is
increased by `options.growth_factor`, up to `options.Δsmax`.

If the number of iterations exceeds the target, the step size is
reduced by `options.shrink_factor`, down to `options.Δsmin`.

Otherwise, the step size remains unchanged.

# Arguments
- `Δs_current::Real`: Current continuation step size.
- `ite::Integer`: Number of Newton iterations performed.
- `options::ContinuationOptions`: Options controlling the adaptive
  step-size strategy.

# Returns
The adapted continuation step size.
"""
function adapt_stepsize(
    Δs_current::U,
    ite::Integer,
    options::ContinuationOptions{U}
) where {U<:Real}

    if ite <= options.target_ite - 1

        Δs_current = min(
            Δs_current * options.growth_factor,
            options.Δsmax
        )

    elseif ite >= options.target_ite + 1

        Δs_current = max(
            Δs_current * options.shrink_factor,
            options.Δsmin
        )

    end

    return Δs_current

end