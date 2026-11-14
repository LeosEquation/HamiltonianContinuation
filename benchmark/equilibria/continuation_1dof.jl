
# ============================================================
# Hamiltonian system
# ============================================================

function hamiltonian_1dof(x, params, t)
    q, p, μ = x

    return p^2 / 2 + q^4 / 4 - q^5 / 5 - μ * q^2 / 2
end

function equations_1dof!(dx, x, params, t)
    q, p, μ = x

    dx[1] = p
    dx[2] = μ * q - q^3 + q^4
    dx[3] = zero(μ)

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

x_ini = [1.0, 0.0, 0.0]

Δs = 0.05
maxsteps = 40

params = nothing

# ============================================================
# Equilibrium continuation benchmark
# ============================================================

trial_1dof = @benchmark EquilibriaContinuation(
    hamiltonian_1dof,
    equations_1dof!,
    boundary_conditions!,
    stopping_condition,
    $x_ini,
    $params,
    $Δs,
    $maxsteps,
) samples=100 evals=1

println("EquilibriaContinuation — 1 DOF")
show(stdout, MIME("text/plain"), trial_1dof)
println()