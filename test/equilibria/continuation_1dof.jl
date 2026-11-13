using Test
using HamiltonianContinuation
using TaylorSeries

set_variables("x", numvars=3)

@testset "EquilibriaContinuation: 1 DOF" begin

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
    # Equilibrium continuation
    # ============================================================

    branch = EquilibriaContinuation(
        hamiltonian_1dof,
        equations_1dof!,
        boundary_conditions!,
        stopping_condition,
        x_ini,
        params,
        Δs,
        maxsteps;
    )

    println("Equilibrium continuation completed.")
    println()

    # ============================================================
    # Results
    # ============================================================

    println("Number of equilibrium points: ", size(branch.x, 1))
    println("Number of limit points: ", length(branch.lp))
    println("Number of branch points: ", length(branch.bp))

    println()
    println("Detected limit points:")

    for (i, lp) in enumerate(branch.lp)
        println("Limit point $i:")
        println("  State: ", lp.x)
        println("  Energy: ", lp.E)
        println("  Residual: ", lp.tol)
    end

    println()
    println("Detected branch points:")

    for (i, bp) in enumerate(branch.bp)
        println("Branch point $i:")
        println("  State: ", bp.x)
        println("  Energy: ", bp.E)
        println("  Residual: ", bp.tol)
    end

    println()

    # ============================================================
    # Equilibrium residuals
    # ============================================================

    residuals = zeros(size(branch.x, 1))

    for i in axes(branch.x, 1)
        dx = similar(branch.x[i, :])

        equations_1dof!(dx, branch.x[i, :], params, 0.0)

        residuals[i] = maximum(abs.(dx[1:2]))
    end

    println("Regular equilibrium residuals:")
    println("  Maximum residual: ", maximum(residuals))
    println("  Minimum residual: ", minimum(residuals))

    # ============================================================
    # Tests
    # ============================================================

    @test size(branch.x, 1) == maxsteps

    # Bifurcation detection
    @test length(branch.lp) == 1
    @test length(branch.bp) == 1

    # Limit point
    @test branch.lp[1].x ≈ [2/3, 0.0, 4/27] atol=1e-8
    @test branch.lp[1].E ≈ -4/405 atol=1e-8
    @test branch.lp[1].tol < 1e-8

    # Branch point
    @test branch.bp[1].x ≈ [0.0, 0.0, 0.0] atol=1e-8
    @test branch.bp[1].E ≈ 0.0 atol=1e-8
    @test branch.bp[1].tol < 1e-8

    # Equilibrium residuals
    @test all(residuals .< 1e-8)

end