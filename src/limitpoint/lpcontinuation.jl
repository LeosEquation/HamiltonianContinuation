

"""
    LPContinuation(f!, g, H, lp_ini::HamiltonianLimitPoint, params, Δs, maxsteps, tol, max_iter)

Performs a numerical continuation of equilibrium solutions starting from a previously 
identified fold (limit point), using a pseudo-arclength predictor-corrector method 
adapted to an extended system formulation.

# Arguments
- `f!`: In-place function representing the system of ODEs (compatible with `DifferentialEquations.jl`) 
and accepting `TaylorN` variables for automatic differentiation.
- `g`: A user-supplied constraint or boundary-checking function `g(q, params, t)` that returns `true` 
if the continuation should be stopped (e.g., variable bounds, conserved quantity violations).
- `H`: A Hamiltonian function `H(x, params, t)` returning the energy at point `x` and time `t`. 
Assumes `t = 0.0` for equilibrium computation.
- `lp_ini::LimitPoint`: A structure representing the initial limit point. Contains the equilibrium 
`x` and an associated null vector of the Jacobian.
- `params`: Parameters to be passed to `f!`, `g`, and `H`.
- `Δs`: Step size for the pseudo-arclength continuation.
- `maxsteps`: Maximum number of continuation steps to perform.
- `tol`: Tolerance for the Newton correction step.
- `max_iter`: Maximum number of Newton iterations per step.

# Returns
- `HamiltonianLimitPointBranch`: A structure containing:
    - `x`: Matrix with each row representing an equilibrium point found along the branch.
    - `E`: Vector containing the Hamiltonian energy evaluated at each equilibrium point.
    - `λ`: Matrix with rows containing the eigenvalue spectrum of the Jacobian at each point.

# Method
The function uses an extended system of dimension `2n - 2` to follow folds (limit points) along a 
bifurcation branch. The continuation direction is computed from the nullspace of the Jacobian of 
this extended system. Each step includes:
- A predictor step along the tangent direction.
- A Newton-based corrector step constrained by pseudo-arclength.
- Jacobian reconstruction via automatic differentiation using `TaylorN`.

The continuation stops if the maximum number of steps is reached, if Newton fails to converge, or 
    if the boundary condition `g` is triggered.

# Notes
- The function assumes that `lp_ini` is a high-accuracy fold point (i.e., the nullspace is 
1-dimensional).
- The eigenvalues stored in `λ` can be used to assess stability along the branch.
- The Hamiltonian energy `E` can be used to classify the branch physically or detect turning 
points in conserved quantities.
"""

function LPContinuation(
    f!,
    lims,
    H,
    lp_ini::HamiltonianLimitPoint{U},
    params,
    Δs_ini::U,
    maxsteps::T;
    atol::U=eps(),
    rtol::U=eps(), ite::T=20,
    reduce_step::U=1.e-1,
    adapt_step::T=3,
) where {U<:Real,T<:Integer}

    ##### Setup: dimensions and auxiliary structures #####

    n = length(lp_ini.x)
    ordtup = [ntuple(k -> count(==(k), (i, j)), TaylorSeries.get_numvars()) for i in 1:n, j in 1:n]

    ##### Output arrays #####

    x = Array{U}(undef, maxsteps, n)
    λ = Array{Complex{U}}(undef, maxsteps, n - 2)
    E = zeros(U, maxsteps)  # <-- new energy array
    tol = Array{U}(undef, maxsteps)

    ##### TaylorN variables for automatic differentiation #####

    xaux = TaylorN.(1:n, order=2)
    dx = zero(xaux)

    Δs = Δs_ini
    adapt_ite = 0

    ##### Evaluation arrays for Jacobian and derivatives #####

    Jeval = zeros(U, n - 2, n)
    dxeval = zeros(U, n - 2)

    ##### Initial extended point (x + null vector) #####

    q0 = zeros(U, 2n - 2)
    q0[1:n] .= lp_ini.x
    q0[n+1:end] .= lp_ini.nullvec
    q1 = copy(q0)

    ##### Extended system: Jacobian, RHS, predictor vector #####

    J = zeros(U, 2n - 2, 2n - 2)
    F = zeros(U, 2n - 2)
    Φ = zeros(U, 2n - 2)
    v = zeros(U, 2n - 2)
    v[end] = one(U)
    Δq = zeros(U, 2n - 2)

    ##### First evaluation step #####

    x[1, :] .= q0[1:n]

    for var in 1:n
        xaux[var][0][1] = q0[var]
    end

    f!(dx, xaux, params, zero(U))

    for row in 1:n-2
        for col in 1:n
            Jeval[row, col] = dx[row][1][col]
        end
        dxeval[row] = dx[row][0][1]
    end

    λ[1, :] .= eigvals(@view(Jeval[:, 1:n-2]))
    E[1] = H(x[1, :], params, zero(U))  # <-- compute energy at first point

    LPJacobian!(J, Jeval, dx, q0, Φ, ordtup, n)
    LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

    tol[1] = norm(F[1:end-1], Inf)

    NS = nullspace(J)

    if size(NS, 2) == 0
        error("The Jacobian at the initial point has no nullspace. Exiting.")
    elseif size(NS, 2) > 1
        println("corank J = $(size(Jeval, 2) - rank(Jeval))")
        error("Jacobian has corank ≥ 2. Exiting.")
    end

    Φ .= NS

    ##### Continuation loop #####

    step = 2
    while step <= maxsteps

        # Predictor step
        for idx in 1:2n-2
            q1[idx] = q0[idx] + Φ[idx] * Δs
        end

        if lims(q1, params, 0.0)
            # @warn("Variable limits exceeded. Exiting.")
            break
        end

        # Evaluate at predicted point
        for var in 1:n
            xaux[var][0][1] = q1[var]
        end

        f!(dx, xaux, params, 0.0)

        for row in 1:n-2
            for col in 1:n
                Jeval[row, col] = dx[row][1][col]
            end
            dxeval[row] = dx[row][0][1]
        end

        # Corrector: Newton iteration
        LPJacobian!(J, Jeval, dx, q1, Φ, ordtup, n)
        LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

        if abs(det(J)) < (atol + rtol * norm(q1, Inf))
            # @warn("System Jacobian was been singular. Stopping continuation.")
            break
        end

        Δq .= J \ F

        iter = 1

        while iter <= ite && norm(Δq, Inf) > (atol + rtol * norm(q1, Inf))

            q1 .-= Δq

            if lims(q1, params, 0.0)
                # @warn("Variable limits exceeded. Exiting.")
                # iter = ite + 1
                break
            end

            for var in 1:n
                xaux[var][0][1] = q1[var]
            end

            f!(dx, xaux, params, 0.0)

            for row in 1:n-2
                for col in 1:n
                    Jeval[row, col] = dx[row][1][col]
                end
                dxeval[row] = dx[row][0][1]
            end

            LPJacobian!(J, Jeval, dx, q1, Φ, ordtup, n)
            LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

            if abs(det(J)) < (atol + rtol * norm(q1, Inf))
                # @warn("System Jacobian was been singular. Stopping continuation.")
                break
            end

            Δq .= J \ F

            iter += 1
        end

        if norm(Δq, Inf) <= (atol + rtol * norm(q1, Inf))
            # Store current solution
            x[step, :] .= q1[1:n]
            λ[step, :] .= eigvals(@view(Jeval[:, 1:n-2]))
            E[step] = H(x[step, :], params, zero(U))  # <-- compute energy at each accepted step
            tol[step] = norm(Δq, Inf)
            # Update direction
            Φ .= J \ v
            normalize!(Φ)
            q0 .= q1
            step += 1
            Δs = Δs_ini
            adapt_ite = 0
        elseif adapt_ite <= adapt_step
            Δs = reduce_step * Δs
            adapt_ite += 1
        else
            @warn("Newton iteration did not converge. Exiting.")
            break
        end
    end

    ##### Return branch of limit points #####

    return HamiltonianLimitPointBranch(x[1:step-1, :], E[1:step-1], λ[1:step-1, :], tol[1:step-1])
end