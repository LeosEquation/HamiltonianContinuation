using HamiltonianContinuation
using TaylorSeries

set_variables("x", numvars = 3)
# ============================================================
# Hamiltonian system
# ============================================================

function pitchfork_hamiltonian(x, params, t)
    q, p, μ = x

    return p^2 / 2 + q^4 / 4 - μ * q^2 / 2
end

# Equilibrium equations (Hamilton's equations)
function pitchfork_equations!(dx, x, params, t)
    q, p, μ = x

    dx[1] = p
    dx[2] = μ * q - q^3
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

x_ini = [0.0, 0.0, -1.0]

Δs = 0.07
maxsteps = 40

params = nothing

# ============================================================
# Equilibrium continuation
# ============================================================

branch = EquilibriaContinuation(
    pitchfork_hamiltonian,
    pitchfork_equations!,
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
println("Detected branch points:")
for (i, bp) in enumerate(branch.bp)
    println("Branch point $i:")
    println("  State: ", bp.x)
    println("  Energy: ", bp.E)
end