function FPContinuation(f!::Function,
    H::Function,
    bc!::Function,
    lims::Function,
    x_ini::Array{U,2},
    v_ini::Array{U,2},
    params,
    Δs_ini::U,
    Δs_lims::Tuple{U,U},
    maxsteps::T,
    Npoints::T,
    integtol::U,
    integorder::T;
    integmaxsteps::T=500,
    parse_eqs::Bool=true,
    atol::U=1.e-12,
    rtol::U=1.e-12,
    ite::T=12,
    nres_max::U=1.e-1) where {U<:Real,T<:Integer}

    println("Running Flip Point Continuation")

    ###

    m1, n = size(x_ini)
    nλ = 3
    nξ = n - nλ
    m2 = Npoints
    m = m1 + m2 * (m1 - 1)
    Nξ = m * nξ
    N = Nξ + nλ

    ###

    Δs_min, Δs_max = Δs_lims

    ###

    ite_1 = ite ÷ 3
    ite_2 = 2 * ite_1

    ###

    nξFP = 2nξ
    NξFP = m * nξFP
    nFP = nξFP + nλ
    NFP = NξFP + nλ

    ###

    progress_step = max(1, floor(Int, maxsteps * 0.001))
    percentage = 0.0

    ###

    Δs = Δs_ini

    ###

    Δτ = one(U) / (m - 1)
    τ = LinRange(zero(U), one(U), m)

    ###

    x = Array{U,3}(undef, maxsteps, m, n)
    μ = Array{Complex{U},2}(undef, maxsteps, nξ)
    v = Array{U,3}(undef, maxsteps, m, nξ)
    E = Array{U,1}(undef, maxsteps)
    tol = Array{U,1}(undef, maxsteps)
    step_size = Array{U,1}(undef, maxsteps)

    ###

    xTFP = TaylorN.([j for i in 1:m-1, j in 1:nFP], order=1)
    Mvec = [zeros(U, nξ, nξ) for i in 1:m-1]


    ###

    x0 = zeros(U, NFP)
    E0 = H(x_ini[1, :], params, 0.0)
    Φ0 = zeros(U, NFP)
    dx0 = zeros(NFP)

    ###

    x_temp = zeros(U, nFP)
    comp2 = 0
    for comp1 in 1:nξ
        comp2 += 1
        x_temp[comp2] = x_ini[1, comp1]
    end
    for comp1 in 1:nξ
        comp2 += 1
        x_temp[comp2] = v_ini[1, comp1]
    end
    for comp1 in 1:nλ
        x_temp[nξFP+comp1] = x_ini[1, nξ+comp1]
    end

    cache_ini = TaylorIntegration.init_cache_optim(zero(U), x_temp, integorder, f!, params; parse_eqs)

    comp3 = 0

    for comp1 in 1:m1-1
        for comp2 in 1:nξ
            comp3 += 1
            x0[comp3] = x_ini[comp1, comp2]
        end
        for comp2 in 1:nξ
            comp3 += 1
            x0[comp3] = v_ini[comp1, comp2]
        end
        # x_temp .= x_ini[comp1, :]
        comp4 = 0
        for comp2 in 1:nξ
            comp4 += 1
            x_temp[comp4] = x_ini[comp1, comp2]
        end
        for comp2 in 1:nξ
            comp4 += 1
            x_temp[comp4] = v_ini[comp1, comp2]
        end
        for comp2 in 1:nλ
            comp4 += 1
            x_temp[comp4] = x_ini[comp1, nξ+comp2]
        end
        for comp2 in 1:m2
            taylorinteg_wrap_optim!(f!, bc!, lims, x_temp, 0.0, Δτ, integtol, cache_ini, params; maxsteps=integmaxsteps)
            for comp4 in 1:nξFP
                comp3 += 1
                x0[comp3] = x_temp[comp4]
            end
        end
    end

    for comp2 in 1:nξ
        comp3 += 1
        x0[comp3] = x_ini[m1, comp2]
    end
    for comp2 in 1:nξ
        comp3 += 1
        x0[comp3] = v_ini[m1, comp2]
    end

    for comp1 in 1:nλ
        x0[NξFP+comp1] = x_ini[m1, nξ+comp1]
    end

    ###

    x1 = copy(x0)
    E1 = copy(E0)
    Φ1 = copy(Φ0)

    ###

    JFP = zeros(U, NFP, NFP)
    FFP = zeros(U, NFP)
    xFP_temp = zeros(U, nFP)
    ΔxFP = zeros(U, NFP)
    vFP = zeros(U, NFP)
    vFP[NFP] = one(U)

    ###

    comp3 = 0

    for comp1 in 1:m
        for comp2 in 1:nξ
            comp3 += 1
            x[1, comp1, comp2] = x0[comp3]
        end
        for comp2 in 1:nξ
            comp3 += 1
            v[1, comp1, comp2] = x0[comp3]
        end
        for comp2 in 1:nλ
            x[1, comp1, nξ+comp2] = x0[NξFP+comp2]
        end
    end

    E[1] = E0

    update_xT!(xTFP, x1, nξFP, nλ, nFP, m, NξFP)

    for comp1 in 1:m

        dx0_view = @view dx0[((comp1-1)*nξFP+1):((comp1-1)*nξFP+nFP)]

        for comp2 in 1:nξ
            xFP_temp[comp2] = x[1, comp1, comp2]
            xFP_temp[nξ+comp2] = v[1, comp1, comp2]
        end
        for comp2 in 1:nλ
            xFP_temp[nξFP+comp2] = x[1, comp1, nξ+comp2]
        end

        f!(dx0_view, xFP_temp, params, τ[comp1])

    end

    ###

    caches = [TaylorIntegration.init_cache_optim(zero(U), xTFP[i, :], integorder, f!, params; parse_eqs) for i in 1:m-1]
    t_test = [true for i in 1:m-1]

    ###

    step = 2

    Threads.@threads for comp1 in 1:m-1
        @views xTview = xTFP[comp1, :]
        t_test[comp1] = taylorinteg_wrap_optim!(f!, bc!, lims, xTview, 0.0, Δτ, integtol, caches[comp1], params; maxsteps=integmaxsteps)
    end

    FPJacobian!(JFP, xTFP, x1, dx0, Φ0, Δτ, nξ, nξFP, nλ, m, NξFP, NFP)
    FPSystem!(FFP, xTFP, x1, x0, dx0, Φ0, Δτ, Δs, nξ, nξFP, NξFP, m, NFP)

    for comp1 in 0:m-1

        FFP_view = @view FFP[comp1*(nξFP)+1:comp1*(nξFP)+nFP]

        bc!(FFP_view, params, τ[comp1+1])

    end

    tol[1] = norm(FFP[1:N-2], Inf)

    NS = nullspace(JFP)

    if size(NS, 2) == 0
        error("Initial point does not have a valid branch direction. Cannot continue.")
    end

    if size(NS, 2) > 1
        @show NS
        error("Initial point is a codimension-2 bifurcation. Cannot continue.")
    end

    Φ0 .= NS

    JFP = nothing

    JFP = spzeros(U, NFP, NFP)

    normΦ = sqrt(Δτ * (0.5 * sum(Φ0[i]^2 for i in 1:nξFP) +
                       sum(Φ0[i]^2 for i in (nξFP+1):(NξFP-nξFP)) +
                       0.5 * sum(Φ0[i]^2 for i in (NξFP-nξFP+1):(NξFP))) +
                 Φ0[NFP-2]^2 + Φ0[NFP-1]^2 + Φ0[NFP]^2)

    Φ0 ./= normΦ

    update_Mvec!(Mvec, xTFP, nξ, m)

    pS = pschur!(Mvec, :L)

    μ[1, :] .= pS.values

    # E1 = H(x1, params, zero(U))

    if sign(Φ0[end]) < 0 && sign(Δs_ini) < 0
        Δs_ini = abs(Δs_ini)
    elseif sign(Φ0[end]) < 0 && sign(Δs_ini) > 0
        Δs_ini = -abs(Δs_ini)
    end

    if sign(Φ0[end]) < 0 && sign(Δs) < 0
        Δs = abs(Δs)
    elseif sign(Φ0[end]) < 0 && sign(Δs) > 0
        Δs = -abs(Δs)
    end

    step_size[1] = abs(Δs)

    real_tol = maximum([atol, rtol * norm(x0, Inf)])

    ###

    while step <= maxsteps

        for comp1 in 1:N
            x1[comp1] = x0[comp1] + Δs * Φ0[comp1]
        end

        update_xT!(xTFP, x1, nξFP, nλ, nFP, m, NξFP)

        Threads.@threads for comp1 in 1:m-1
            @views xTview = xTFP[comp1, :]
            t_test[comp1] = taylorinteg_wrap_optim!(f!, bc!, lims, xTview, 0.0, Δτ, integtol, caches[comp1], params; maxsteps=integmaxsteps)
        end

        FPJacobian!(JFP, xTFP, x1, dx0, Φ0, Δτ, nξ, nξFP, nλ, m, NξFP, NFP)
        FPSystem!(FFP, xTFP, x1, x0, dx0, Φ0, Δτ, Δs, nξ, nξFP, NξFP, m, NFP)

        for comp1 in 0:m-1

            FFP_view = @view FFP[comp1*(nξFP)+1:comp1*(nξFP)+nFP]

            bc!(FFP_view, params, τ[comp1+1])

        end

        nres = norm(FFP, Inf)
        nstep = Inf # norm(ΔxFP, Inf)

        newton_iter = 1

        while all(t_test) &&
                  newton_iter <= ite &&
                  nres <= nres_max &&
                  (nres > real_tol ||
                   nstep > real_tol)

            ΔxFP .= JFP \ FFP

            x1 .-= ΔxFP

            update_xT!(xTFP, x1, nξFP, nλ, nFP, m, NξFP)

            Threads.@threads for comp1 in 1:m-1
                @views xTview = xTFP[comp1, :]
                t_test[comp1] = taylorinteg_wrap_optim!(f!, bc!, lims, xTview, 0.0, Δτ, integtol, caches[comp1], params; maxsteps=integmaxsteps, reltol=0.0)
            end

            FPJacobian!(JFP, xTFP, x1, dx0, Φ0, Δτ, nξ, nξFP, nλ, m, NξFP, NFP)
            FPSystem!(FFP, xTFP, x1, x0, dx0, Φ0, Δτ, Δs, nξ, nξFP, NξFP, m, NFP)

            for comp1 in 0:m-1

                FFP_view = @view FFP[comp1*(nξFP)+1:comp1*(nξFP)+nFP]

                bc!(FFP_view, params, τ[comp1+1])

            end

            nres = norm(FFP, Inf)
            nstep = norm(ΔxFP, Inf)

            newton_iter += 1

        end

        if all(t_test) && (nres <= real_tol || nstep <= real_tol)

            tol[step] = nstep

            step_size[step] = abs(Δs)

            update_Mvec!(Mvec, xTFP, nξ, m)

            pS = pschur!(Mvec, :L)

            μ[step, :] .= pS.values

            comp3 = 0

            for comp1 in 1:m
                for comp2 in 1:nξ
                    comp3 += 1
                    x[step, comp1, comp2] = x1[comp3]
                end
                for comp2 in 1:nξ
                    comp3 += 1
                    v[step, comp1, comp2] = x1[comp3]
                end
                for comp2 in 1:nλ
                    x[step, comp1, nξ+comp2] = x1[NξFP+comp2]
                end
            end

            E1 = H(x[step, 1, :], params, τ[1])
            E[step] = E1

            for comp1 in 1:m

                dx0_view = @view dx0[((comp1-1)*nξFP+1):((comp1-1)*nξFP+nFP)]

                for comp2 in 1:nξ
                    xFP_temp[comp2] = x[step, comp1, comp2]
                    xFP_temp[nξ+comp2] = v[step, comp1, comp2]
                end
                for comp2 in 1:nλ
                    xFP_temp[nξFP+comp2] = x[step, comp1, nξ+comp2]
                end

                f!(dx0_view, xFP_temp, params, τ[comp1])

            end

            Φ1 .= JFP \ vFP

            normΦ = sqrt(Δτ * (0.5 * sum(Φ1[i]^2 for i in 1:nξFP) +
                               sum(Φ1[i]^2 for i in (nξFP+1):(NξFP-nξFP)) +
                               0.5 * sum(Φ1[i]^2 for i in (NξFP-nξFP+1):(NξFP))) +
                         Φ1[NFP-2]^2 + Φ1[NFP-1]^2 + Φ1[NFP]^2)

            Φ1 ./= normΦ

            x0 .= x1
            Φ0 .= Φ1
            E0 = E1
            real_tol = maximum([atol, rtol * norm(x0, Inf)])

            if newton_iter <= ite_1
                Δs *= 1.5
            elseif ite_2 ≤ newton_iter < ite
                Δs *= 0.5
            end

            Δs = sign(Δs) * clamp(abs(Δs), Δs_min, Δs_max)

            if step % progress_step == 0 || step == maxsteps
                percentage = round(100.0 * step / maxsteps, digits=2)

                print(" \r Progress : $percentage % \t                       ")
            end

            step += 1

        else

            Δs *= 0.25

            if abs(Δs) ≤ Δs_min

                println("\n-------------------------")
                @warn("System norm did not converge to tolerance. Stopping continuation.")
                println("-------------------------")
                break

            end

        end

    end

    println("")

    return HamiltonianFlipPointBranch(τ,
        x[1:step-1, :, :],
        v[1:step-1, :, :],
        μ[1:step-1, :],
        E[1:step-1],
        tol[1:step-1],
        step_size[1:step-1])

end