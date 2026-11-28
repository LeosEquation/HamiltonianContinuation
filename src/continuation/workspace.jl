"""
    ContinuationWorkspace{B,S,T,P,F}

Stores the data structures required during a continuation process.

# Fields
- `branch::B`: Stores the continuation branch.
- `state::S`: Stores the current numerical state and auxiliary data required
  by the continuation algorithm.
- `taylor::T`: Stores the Taylor expansions of the position and velocity
  variables used during the continuation process.
- `phase_space::P`: Stores the position and velocity variables of the
  physical system, independently of their representation.
- `bifurcations::F`: Stores the bifurcation data associated with the
  continuation process.

The concrete types and internal representations of these fields depend
on the continuation method.
"""
struct ContinuationWorkspace{B,S,T,P,F}
    branch::B
    state::S
    taylor::T
    phase_space::P
    bifurcations::F
end

"""
    build_continuation_workspace(::Type{C}, problem, options)

Builds the workspace required to perform a continuation process.

The workspace components are constructed using `build_continuation_state`,
`build_branch`, `build_taylor_workspace`, `build_phase_space`, and
`build_bifurcations`, which dispatch according to the continuation type
`C` and the problem.

# Arguments
- `C::Type{<:AbstractContinuation}`: Continuation method to be used.
- `problem::AbstractProblem`: Problem defining the system to be continued.
- `options::ContinuationOptions`: Continuation algorithm options.

# Returns
A tuple containing:
- `workspace::ContinuationWorkspace`: Workspace containing the continuation
  branch, numerical state, Taylor expansions, phase-space variables, and
  bifurcation data.
- `Δs_current`: Initial continuation step size, initialized from `options.Δs`.
- `maxsteps`: Maximum number of continuation steps, obtained from
  `options.maxsteps`.
"""
function build_continuation_workspace(
    ::Type{C},
    problem::P,
    options::ContinuationOptions{U}
) where {C<:AbstractContinuation,P<:AbstractProblem,U<:Real}

    (; Δs, maxsteps) = options

    taylor = build_taylor_workspace(problem, space)

    phase_space = build_phase_space(C, taylor, problem)

    state = build_continuation_state(C, taylor, problem)

    branch = build_branch(C, phase_space, problem, options)

    bifurcations = build_bifurcations(C, problem)

    Δs_current = copy(Δs)

    workspace = ContinuationWorkspace(
        branch,
        state,
        taylor,
        phase_space,
        bifurcations
    )

    return workspace, Δs_current, maxsteps
end