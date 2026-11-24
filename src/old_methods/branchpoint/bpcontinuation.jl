
"""
    BPContinuation(f!, g, H, bp_ini::BranchPoint, params, Δs, maxsteps, tol, max_iter)

This version also evaluates the Hamiltonian energy at each equilibrium point.

# Arguments
- `H`: Hamiltonian function `H(x, params, t)` that returns the energy at a given point `x` 
and time `t`. Since all equilibria are time-independent, `t = 0.0` is assumed for all evaluations.

# Returns
- `HamiltonianBranchPointBranch`: A structure containing:
    - `x`: Matrix where each row is an equilibrium point found during continuation.
    - `E`: Vector with Hamiltonian energy evaluated at each equilibrium.
    - `λ`: Matrix containing the eigenvalue spectra (from the reduced Jacobian) at each point.
    - `β`: Placeholder vector (initialized as zeros) for future classification of the branch points.

# Notes
- Energy values in `E` allow for global energetic interpretation of the solution branch.
"""
function BPContinuation(
    f!,
    lims,
    H,
    bp_ini::HamiltonianBranchPoint{U},
    params,
    Δs_ini::U,
    maxsteps::T;
    atol::U=eps(),
    rtol::U=eps(),
    ite::T=20,
    reduce_step::U=1.e-1,
    adapt_step::T=3,
) where {U<:Real,T<:Integer}

    # --- Setup ---
    n = length(bp_ini.x)
    ordtup = [ntuple(k -> count(==(k), (i, j)), TaylorSeries.get_numvars()) for i in 1:n, j in 1:n]

    x = Array{U}(undef, maxsteps, n)
    λ = Array{Complex{U}}(undef, maxsteps, n - 2)
    E = Array{U}(undef, maxsteps)
    β = zeros(U, maxsteps)
    tol = Array{U}(undef, maxsteps)

    Δs = Δs_ini
    adapt_ite = 0

    δx = TaylorN.(1:n, order=2)
    xaux = copy(δx)
    dx = zero(δx)

    # Initial condition q0
    q0 = zeros(U, 2n - 1)
    q0[1:n] .= bp_ini.x
    q0[n+1:2n-2] .= bp_ini.nullvec
    q1 = copy(q0)

    # Jacobian system variables
    J = zeros(U, 2n - 1, 2n - 1)
    F = zeros(U, 2n - 1)
    Φ = zeros(U, 2n - 1)
    v = zeros(U, 2n - 1)
    v[end] = one(U)
    Δq = zeros(U, 2n - 1)

    # Reduced Jacobian
    Jeval = zeros(U, n - 2, n)
    dxeval = zeros(U, n - 2)

    # Store first point
    x[1, :] .= q0[1:n]
    E[1] = H(@view(q0[1:n]), params, zero(U))

    for j in 1:n
        xaux[j][0][1] = q0[j]
    end

    f!(dx, xaux, params, zero(U))

    for row in 1:n-2
        for col in 1:n
            Jeval[row, col] = dx[row][1][col]
        end
        dxeval[row] = dx[row][0][1]
    end

    BPJacobian!(J, Jeval, q0, dx, Φ, ordtup, n)
    BPSystem!(F, dxeval, Jeval, q1, q0, Φ, Δs, n)

    tol[1] = norm(F[1:end-1], Inf)

    NS = nullspace(J)
    if size(NS, 2) == 0
        error("No nullspace direction found. Aborting.")
    elseif size(NS, 2) > 1
        error("More than one nullspace direction found. Aborting.")
    end
    Φ .= NS
    λ[1, :] .= eigvals(@view(Jeval[:, 1:n-2]))

    step = 2
    while step <= maxsteps

        for idx in 1:2n-1
            q1[idx] = q0[idx] + Φ[idx] * Δs
        end

        if lims(q1, params, 0.0)
            # @warn("Initial System Jacobian was been singular. Aborting.")
            break
        end

        for j in 1:n
            xaux[j][0][1] = q1[j]
        end

        f!(dx, xaux, params, zero(U))

        for row in 1:n-2
            for col in 1:n
                Jeval[row, col] = dx[row][1][col]
            end
            dxeval[row] = dx[row][0][1]
        end

        BPJacobian!(J, Jeval, q1, dx, Φ, ordtup, n)
        BPSystem!(F, dxeval, Jeval, q1, q0, Φ, Δs, n)

        if abs(det(J)) < (atol + rtol * norm(q1, Inf))
            # @warn("Initial System Jacobian was been singular. Aborting.")
            break
        end

        Δq .= J \ F

        iter = 1

        while iter <= ite && norm(Δq, Inf) > (atol + rtol * norm(q1, Inf))

            q1 .-= Δq

            if lims(q1, params, 0.0)
                # @warn "Constraint `g` triggered during Newton. Exiting."
                break
            end

            for j in 1:n
                xaux[j][0][1] = q1[j]
            end

            f!(dx, xaux, params, zero(U))

            for row in 1:n-2
                for col in 1:n
                    Jeval[row, col] = dx[row][1][col]
                end
                dxeval[row] = dx[row][0][1]
            end

            BPJacobian!(J, Jeval, q1, dx, Φ, ordtup, n)
            BPSystem!(F, dxeval, Jeval, q1, q0, Φ, Δs, n)

            if abs(det(J)) < (atol + rtol * norm(q1, Inf))
                # @warn("System Jacobian was been singular. Stopping continuation.")
                break
            end

            Δq .= J \ F

            iter += 1
        end

        if norm(Δq, Inf) <= (atol + rtol * norm(q1, Inf))
            x[step, :] .= q1[1:n]
            λ[step, :] .= eigvals(@view(Jeval[:, 1:n-2]))
            E[step] = H(@view(q1[1:n]), params, zero(U))
            tol[step] = norm(Δq, Inf)
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

    return HamiltonianBranchPointBranch(x[1:step-1, :], E[1:step-1], λ[1:step-1, :], β[1:step-1], tol[1:step-1])
end
