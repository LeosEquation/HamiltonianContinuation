@testset "EquilibriaContinuation: 2 DOF" begin

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
    # Equilibrium continuation
    # ============================================================

    branch = EquilibriaContinuation(
        hamiltonian_2dof,
        equations_2dof!,
        boundary_conditions!,
        stopping_condition,
        x_ini,
        params,
        Δs,
        maxsteps;
    )

    # ============================================================
    # Equilibrium residuals
    # ============================================================

    residuals = zeros(size(branch.x, 1))

    for i in axes(branch.x, 1)
        dx = similar(branch.x[i, :])

        equations_2dof!(dx, branch.x[i, :], params, 0.0)

        residuals[i] = maximum(abs.(dx[1:4]))
    end

    # ============================================================
    # Tests
    # ============================================================

    @test size(branch.x, 1) == maxsteps

    # Bifurcation detection
    @test length(branch.lp) == 1
    @test length(branch.bp) == 1

    # Limit point
    @test branch.lp[1].x ≈ [2/3, 0.0, 0.0, 0.0, 4/27] atol=1e-8
    @test branch.lp[1].E ≈ -4/405 atol=1e-8
    @test branch.lp[1].tol < 1e-8

    # Branch point
    @test branch.bp[1].x ≈ [0.0, 0.0, 0.0, 0.0, 0.0] atol=1e-8
    @test branch.bp[1].E ≈ 0.0 atol=1e-8
    @test branch.bp[1].tol < 1e-8

    # Equilibrium residuals
    @test all(residuals .< 1e-8)

end