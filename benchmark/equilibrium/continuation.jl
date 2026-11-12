
using BenchmarkTools
using HamiltonianContinuation
using TaylorSeries

set_variables("x", numvars=3)

# ============================================================
# Pitchfork system
# ============================================================

function pitchfork_hamiltonian(x, params, t)
    q, p, μ = x

    return p^2 / 2 + q^4 / 4 - μ * q^2 / 2
end

function pitchfork_equations!(dx, x, params, t)
    q, p, μ = x

    dx[1] = p
    dx[2] = μ * q - q^3
    dx[3] = zero(μ)

    return nothing
end

# ============================================================
# Limit point system
# ============================================================

function limitpoint_hamiltonian(x, params, t)
    q, p, μ = x

    return p^2 / 2 + q^3 / 3 - μ * q
end

function limitpoint_equations!(dx, x, params, t)
    q, p, μ = x

    dx[1] = p
    dx[2] = μ - q^2
    dx[3] = zero(μ)

    return nothing
end

# ============================================================
# Common functions
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

params = nothing

maxsteps = 40

x_ini_pitchfork = [0.0, 0.0, -1.0]
Δs_pitchfork = 0.07

x_ini_limitpoint = [1.0, 0.0, 1.0]
Δs_limitpoint = -0.07

# ============================================================
# Performance benchmarks
# ============================================================

println("Benchmarking EquilibriaContinuation...")
println()

# Pitchfork (BP)
println("Pitchfork (BP):")

pitchfork_trial = @benchmark EquilibriaContinuation(
    $pitchfork_hamiltonian,
    $pitchfork_equations!,
    $boundary_conditions!,
    $stopping_condition,
    $x_ini_pitchfork,
    $params,
    $Δs_pitchfork,
    $maxsteps
)

display(pitchfork_trial)

println()
println("========================================")
println()

# Limit point (LP)
println("Limit point (LP):")

limitpoint_trial = @benchmark EquilibriaContinuation(
    $limitpoint_hamiltonian,
    $limitpoint_equations!,
    $boundary_conditions!,
    $stopping_condition,
    $x_ini_limitpoint,
    $params,
    $Δs_limitpoint,
    $maxsteps
)

display(limitpoint_trial)

println()
println("Benchmark completed.")
