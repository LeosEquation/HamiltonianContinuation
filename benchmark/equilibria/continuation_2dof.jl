# ============================================================
# Hamiltonian system
# ============================================================

function hamiltonian_2dof(x, params, t)
    q1, p1, q2, p2, μ = x

    return p1^2 / 2 + q1^4 / 4 - q1^5 / 5 - μ * q1^2 / 2 +
        p2^2 / 2 + q2^2 / 2
end

function equations_2dof!(dx, x, params, t)
    q1, p1, q2, p2, μ = x

    dx[1] = p1
    dx[2] = μ * q1 - q1^3 + q1^4
    dx[3] = p2
    dx[4] = -q2
    dx[5] = zero(μ)

    return nothing
end

# ============================================================
# Boundary and stopping conditions
# ============================================================

function boundary_conditions!(x, params, t)
    return nothing
end

function stopping_condition(x, params, t)
    return false
end

# ============================================================
# Initial conditions and continuation parameters
# ============================================================

x_ini = [1.0, 0.0, 0.0, 0.0, 0.0]

Δs = 0.05
maxsteps = 40

params = nothing

# ============================================================
# Equilibrium continuation benchmark
# ============================================================

trial_2dof = @benchmark EquilibriaContinuation(
    hamiltonian_2dof,
    equations_2dof!,
    boundary_conditions!,
    stopping_condition,
    $x_ini,
    $params,
    $Δs,
    $maxsteps,
) samples=100 evals=1

println("EquilibriaContinuation — 2 DOF")
show(stdout, MIME("text/plain"), trial_2dof)
println()
