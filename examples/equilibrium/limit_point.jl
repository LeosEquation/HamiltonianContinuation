using Revise
using HamiltonianContinuation
using TaylorSeries

set_variables("x", numvars=3)

# ============================================================
# Hamiltonian system
# ============================================================

function limitpoint_hamiltonian(x, params, t)
    q, p, μ = x

    return p^2 / 2 + q^3 / 3 - μ * q
end

# Equilibrium equations (Hamilton's equations)
function limitpoint_equations!(dx, x, params, t)
    q, p, μ = x

    dx[1] = p
    dx[2] = μ - q^2
    dx[3] = zero(μ)

    return nothing
end

# No additional boundary conditions
function boundary_conditions!(x, params, t)
    return nothing
end

# No stopping condition
function stopping_condition(x, params, t)
    return false
end

# ============================================================
# Initial conditions and continuation parameters
# ============================================================

x_ini = [1.0, 0.0, 1.0]

Δs = -0.07
maxsteps = 40

params = nothing

# ============================================================
# Equilibrium continuation
# ============================================================

branch = EquilibriaContinuation(
    limitpoint_hamiltonian,
    limitpoint_equations!,
    boundary_conditions!,
    stopping_condition,
    x_ini,
    params,
    Δs,
    maxsteps;
)

# ============================================================
# Results
# ============================================================

println("Equilibrium continuation completed.")
println()

println("Number of equilibrium points: ", size(branch.x, 1))
println("Number of limit points: ", length(branch.lp))
println("Number of branch points: ", length(branch.bp))

println()
println("Final equilibrium point:")
println(branch.x[end, :])

println()
println("Detected limit points:")
for (i, lp) in enumerate(branch.lp)
    println("Limit point $i:")
    println("  State: ", lp.x)
    println("  Energy: ", lp.E)
end