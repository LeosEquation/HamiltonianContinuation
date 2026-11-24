function OrbitContinuation(
    f!,
    g!,
    H,
    bc!,
    lims,
    x_ini::Array{U,2},
    params,
    Δs_ini::U,
    maxsteps::T,
    Npoints::T,
    integtol::U,
    integorder::T;
    integmaxsteps::T=500,
    parse_eqs::Bool=true,
    atol::U=eps(),
    rtol::U=eps(),
    ite::T=20,
    aite::T=3,) where {U<:Real,T<:Integer}

    println("Running Orbit Continuation.")

    ###

    m1, n = size(x_ini)
    nξ = n - 2
    nλ = 2
    m2 = Npoints
    m = m1 + m2 * (m1 - 1)
    Nξ = m * nξ
    N = Nξ + nλ

    progress_step = max(1, floor(Int, maxsteps * 0.001))
    percentage = 0.0

    Δτ = one(U) / (m - 1)
    τ = LinRange(zero(U), one(U), m)

    ###

    Δs = Δs_ini
    astep = 0

    ###

    x = Array{U,3}(undef, maxsteps, m, n)
    μ = Array{Complex{U},2}(undef, maxsteps, n - 2)
    E = Array{U,1}(undef, maxsteps)
    stbl = Array{Bool,1}(undef, maxsteps)
    tol = Array{U,1}(undef, maxsteps)
    FP = HamiltonianFlipPoint{U}[]
    BCP = HamiltonianBranchCyclePoint{U}[]
    FPtest = Array{U,1}(undef, maxsteps)
    BCPtest = Array{U,1}(undef, maxsteps)

    ###

    xT = TaylorN.([j for i in 1:m-1, j in 1:n], order=1)
    Mvec = [zeros(U, n - 2, n - 2) for i in 1:m-1]

    ###

    xHTN = TaylorN.(1:n, order=1)
    HTN = TaylorN(1, order=1)

    ###

    x0 = zeros(U, N)
    E0 = 0.0
    FPtest0 = 0.0
    BCPtest0 = 0.0

    ###

    x_temp = zeros(U, n)
    for comp1 in 1:nξ
        x_temp[comp1] = x_ini[1, comp1]
    end
    for comp1 in 1:nλ
        x_temp[nξ+comp1] = x_ini[1, nξ+comp1]
    end

    cache_ini = TaylorIntegration.init_cache_optim(zero(U), x_temp, integorder, f!, params; parse_eqs)

    comp3 = 0

    for comp1 in 1:m1-1
        for comp2 in 1:nξ
            comp3 += 1
            x0[comp3] = x_ini[comp1, comp2]
        end
        x_temp .= x_ini[comp1, :]
        for comp2 in 1:m2
            taylorinteg_wrap_optim!(f!, bc!, lims, x_temp, 0.0, Δτ, integtol, cache_ini, params; maxsteps=integmaxsteps)
            for comp2 in 1:nξ
                comp3 += 1
                x0[comp3] = x_temp[comp2]
            end
        end
    end

    for comp2 in 1:nξ
        comp3 += 1
        x0[comp3] = x_ini[m1, comp2]
    end

    for comp1 in 1:nλ
        x0[Nξ+comp1] = x_ini[1, nξ+comp1]
    end

    E0 = H(x_ini[1, :], params, 0.0)

    E[1] = E0

    Φ0 = zeros(U, N)

    dx0 = zeros(N)

    ###

    x1 = copy(x0)
    E1 = 1.0
    FPtest1 = 0.0
    BCPtest1 = 0.0
    Φ1 = copy(Φ0)

    ###

    J = zeros(U, N, N)
    F = zeros(U, N)
    Δx = zeros(U, N)
    v = zeros(U, N)
    v[N] = one(U)

    ###

    nFP = 2 * nξ + nλ
    nξFP = 2 * nξ
    NFP = m * nξFP + nλ
    NξFP = m * nξFP
    xTFP = TaylorN.([j for i in 1:m-1, j in 1:nFP], order=1)
    xFP = zeros(U, NFP)
    x0FP = zeros(U, NFP)
    ΔxFP = zeros(U, NFP)
    xFP_temp = zeros(U, nFP)
    xFP_eval = zeros(U, m, nFP)
    ξFP_eval = zeros(U, m, n)
    vFP_eval = zeros(U, m, nξ)
    EFP = zero(U)
    dx0FP = zeros(U, NFP)

    ###

    comp4 = 0
    comp3 = 0
    for comp1 in 1:m
        for comp2 in 1:nξ
            comp4 += 1
            comp3 += 1
            xFP[comp4] = x0[comp3]
        end
        for comp2 in 1:nξ
            comp4 += 1
            xFP[comp4] = 1.0
        end
    end

    for comp1 in 1:nλ
        xFP[NξFP+comp1] = x_ini[1, nξ+comp1]
    end

    ### 

    JFP = spzeros(U, NFP, NFP)
    FFP = zeros(U, NFP)

    ### 

    comp2 = 0
    for comp1 in 1:nξ
        comp2 += 1
        xFP_temp[comp2] = x_ini[1, comp1]
    end
    for comp1 in 1:nξ
        comp2 += 1
        xFP_temp[comp2] = 1.0
    end
    for comp1 in 1:nλ
        comp2 += 1
        xFP_temp[comp2] = x_ini[1, nξ+comp1]
    end

    cache_ini_FP = TaylorIntegration.init_cache_optim(zero(U), xFP_temp, integorder, g!, params; parse_eqs)

    ###

    comp3 = 0
    for comp1 in 1:m
        for comp2 in 1:nξ
            comp3 += 1
            x[1, comp1, comp2] = x0[comp3]
        end
        for comp2 in 1:nλ
            x[1, comp1, nξ+comp2] = x0[Nξ+comp2]
        end
    end

    for comp1 in 1:m

        dx0_view = @view dx0[((comp1-1)*nξ+1):((comp1-1)*nξ+n)]

        f!(dx0_view, x[1, comp1, :], params, τ[comp1])

    end

    update_xT!(xT, x0, nξ, nλ, n, m, Nξ)
    update_xT!(xTFP, xFP, nξFP, nλ, nFP, m, NξFP)

    ###

    caches = [TaylorIntegration.init_cache_optim(zero(U), xT[i, :], integorder, f!, params; parse_eqs) for i in 1:m-1]
    caches_FP = [TaylorIntegration.init_cache_optim(zero(U), xTFP[i, :], integorder, g!, params; parse_eqs) for i in 1:m-1]

    ###

    for comp1 in 1:m-1
        @views xTview = xT[comp1, :]
        taylorinteg_wrap_optim!(f!, bc!, lims, xTview, 0.0, Δτ,
            integtol, caches[comp1], params; maxsteps=integmaxsteps)
    end

    OrbitJacobian!(J, xT, dx0, Φ0, Δτ, nξ, nλ, m, Nξ, N)
    OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, Δτ, nξ, m, Nξ, N)

    for comp1 in 1:m

        F_view = @view F[(comp1-1)*nξ+1:(comp1-1)*nξ+n]

        bc!(F_view, params, τ[comp1])

    end

    tol[1] = norm(F[1:N-2], Inf)

    NS = nullspace(J)

    if size(NS, 2) == 0
        error("Initial point does not have a valid branch direction. Cannot continue.")
    end

    if size(NS, 2) > 1
        @show NS
        error("Initial point is a codimension-2 bifurcation. Cannot continue.")
    end

    Φ0 .= NS

    J = nothing

    J = spzeros(U, N, N)

    normΦ = sqrt(Δτ * (0.5 * sum(Φ0[i]^2 for i in 1:nξ) +
                       sum(Φ0[i]^2 for i in (nξ+1):(Nξ-nξ)) +
                       0.5 * sum(Φ0[i]^2 for i in (Nξ-nξ+1):(Nξ))) +
                 Φ0[N-1]^2 + Φ0[N]^2)

    Φ0 ./= normΦ

    update_Mvec!(Mvec, xT, nξ, m)

    pS = pschur!(Mvec, :L)

    BCPtest0 = prod(real(pS.values) .- 1.0)
    BCPtest[1] = BCPtest0

    FPtest0 = prod(real(pS.values) .+ 1.0)
    FPtest[1] = FPtest0

    μ[1, :] .= pS.values

    μ_test = maximum(abs.(μ[1, :]))

    if (μ_test - 1.0) > 0.01

        stbl[1] = false

    else

        stbl[1] = true

    end

    for comp1 in 1:N
        x1[comp1] = x0[comp1] + Δs * Φ0[comp1]
    end

    for comp1 in 1:n
        xHTN[comp1][0][1] = x[1, 1, comp1]
    end

    HTN .= H(xHTN, params, zero(U))

    ΔH = sum(HTN[1][i] * Φ0[i] for i in 1:nξ) + sum(HTN[1][i] * Φ0[Nξ+i] for i in 1:nλ)

    if sign(ΔH) < 0 && sign(Δs_ini) < 0
        Δs_ini = abs(Δs_ini)
    elseif sign(ΔH) < 0 && sign(Δs_ini) > 0
        Δs_ini = -abs(Δs_ini)
    end

    if sign(ΔH) < 0 && sign(Δs) < 0
        Δs = abs(Δs)
    elseif sign(ΔH) < 0 && sign(Δs) > 0
        Δs = -abs(Δs)
    end

    ###

    step = 2

    while step <= maxsteps

        for comp1 in 1:N
            x1[comp1] = x0[comp1] + Δs * Φ0[comp1]
        end

        update_xT!(xT, x1, nξ, nλ, n, m, Nξ)

        for comp1 in 1:m-1
            @views xTview = xT[comp1, :]
            taylorinteg_wrap_optim!(f!, bc!, lims, xTview, 0.0, Δτ,
                integtol, caches[comp1], params; maxsteps=integmaxsteps)
        end

        OrbitJacobian!(J, xT, dx0, Φ0, Δτ, nξ, nλ, m, Nξ, N)
        OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, Δτ, nξ, m, Nξ, N)

        for comp1 in 1:m

            F_view = @view F[(comp1-1)*nξ+1:(comp1-1)*nξ+n]

            bc!(F_view, params, τ[comp1])

        end

        if isnan(det(J)) || abs(det(J)) < (atol + rtol * norm(x1, Inf))
            # error_massage = "System Jacobian was been singular. Stopping continuation."
            break
        end

        Δx .= J \ F

        newton_iter = 1

        while newton_iter <= ite && norm(Δx, Inf) > (atol + rtol * norm(x1, Inf))

            # @show (j, norm(F))

            x1 .-= Δx

            update_xT!(xT, x1, nξ, nλ, n, m, Nξ)

            for comp1 in 1:m-1
                @views xTview = xT[comp1, :]
                taylorinteg_wrap_optim!(f!, bc!, lims, xTview, 0.0, Δτ,
                    integtol, caches[comp1], params; maxsteps=integmaxsteps)
            end

            OrbitJacobian!(J, xT, dx0, Φ0, Δτ, nξ, nλ, m, Nξ, N)
            OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, Δτ, nξ, m, Nξ, N)

            for comp1 in 1:m

                F_view = @view F[(comp1-1)*(n-2)+1:(comp1-1)*(n-2)+n]
                bc!(F_view, params, τ[comp1])

            end

            if isnan(det(J)) || abs(det(J)) < (atol + rtol * norm(x1, Inf))
                # error_massage = "System Jacobian was been singular. Stopping continuation."
                break
            end

            Δx .= J \ F

            newton_iter += 1

        end

        # if any(isnan.(Δx)) || norm(Δx, Inf) > (atol + rtol * norm(x1, Inf))
        #     println("")
        #     @warn("System norm did not converge to tolerance. Stopping continuation.")
        #     break
        # end

        if all(.!isnan.(Δx)) && norm(Δx, Inf) <= (atol + rtol * norm(x1, Inf))

            tol[step] = norm(Δx, Inf)

            comp3 = 0

            for comp1 in 1:m
                for comp2 in 1:nξ
                    comp3 += 1
                    x[step, comp1, comp2] = x1[comp3]
                end
                x[step, comp1, n-1] = x1[N-1]
                x[step, comp1, n] = x1[N]
            end

            for comp1 in 1:m

                dx0_view = @view dx0[((comp1-1)*nξ+1):((comp1-1)*nξ+n)]

                f!(dx0_view, x[step, comp1, :], params, τ[comp1])

            end

            E1 = H(x[step, 1, :], params, τ[1])
            E[step] = E1

            update_Mvec!(Mvec, xT, nξ, m)

            pS = pschur!(Mvec, :L)

            BCPtest1 = real(prod(pS.values .- 1.0))
            BCPtest[step] = BCPtest1

            FPtest1 = real(prod(pS.values .+ 1.0))
            FPtest[step] = FPtest1

            μ[step, :] .= pS.values

            μ_test = maximum(abs.(μ[step, :]))

            if (μ_test - 1.0) > 0.01

                stbl[step] = false

            else

                stbl[step] = true

            end

            Φ1 .= J \ v

            normΦ = sqrt(Δτ * (0.5 * sum(Φ1[j]^2 for j in 1:nξ) +
                               sum(Φ1[j]^2 for j in (nξ+1):(Nξ-nξ)) +
                               0.5 * sum(Φ1[j]^2 for j in (Nξ-nξ+1):(Nξ))) +
                         Φ1[N-1]^2 + Φ1[N]^2)

            Φ1 ./= normΦ

            if sign(FPtest1) != sign(FPtest0)

                print(" \r Progress : $percentage % \t | #FP =... \t |                       ")

                update_Mvec!(Mvec, xT, nξ, m)
                Mtotal = reduce((A, B) -> B * A, Mvec)
                eigM = eigen(Mtotal)
                idx = argmin(abs.(eigM.values .- (-1.0 + 0.0im)))

                # @show E0
                # @show E1
                # @show eigM.values

                comp2 = 0
                for comp1 in 1:nξ
                    comp2 += 1
                    xFP_temp[comp2] = x1[comp1]
                end
                for comp1 in 1:nξ
                    comp2 += 1
                    xFP_temp[comp2] = real(eigM.vectors[comp1, idx])
                end
                for comp1 in 1:nλ
                    comp2 += 1
                    xFP_temp[comp2] = x1[Nξ+comp1]
                end

                comp3 = 0
                for comp1 in 1:m
                    for comp2 in 1:nξFP
                        comp3 += 1
                        xFP[comp3] = xFP_temp[comp2]
                    end
                end

                comp3 = nξFP

                for comp1 in 2:m

                    taylorinteg_wrap_optim!(g!, bc!, lims, xFP_temp, 0.0, Δτ, integtol, cache_ini_FP, params; maxsteps=integmaxsteps)

                    for comp2 in 1:nξFP

                        comp3 += 1

                        xFP[comp3] = xFP_temp[comp2]

                    end

                end

                for comp1 in 1:nλ

                    comp3 += 1

                    xFP[comp3] = xFP_temp[nξFP+comp1]

                end

                x0FP .= xFP

                comp3 = 0
                for comp1 in 1:m
                    for comp2 in 1:nξFP
                        comp3 += 1
                        xFP_eval[comp1, comp2] = xFP[comp3]
                    end
                    for comp2 in 1:nλ
                        xFP_eval[comp1, nξFP+comp2] = xFP[NξFP+comp2]
                    end
                end

                for comp1 in 1:m

                    @views dx0_view = dx0FP[((comp1-1)*nξFP+1):((comp1-1)*nξFP+nFP)]
                    g!(dx0_view, xFP_eval[comp1, :], params, τ[comp1])

                end

                update_xT!(xTFP, xFP, nξFP, nλ, nFP, m, NξFP)

                for comp1 in 1:m-1
                    @views xTview = xTFP[comp1, :]
                    taylorinteg_wrap_optim!(g!, bc!, lims, xTview, 0.0, Δτ,
                        integtol, caches_FP[comp1], params; maxsteps=integmaxsteps)
                end

                FPJacobian!(JFP, xTFP, xFP, dx0FP, Δτ, nξ, nξFP, nλ, m, NξFP, NFP)
                FPSystem!(FFP, xTFP, xFP, x0FP, dx0FP, Δτ, nξ, nξFP, m, NFP)

                for comp1 in 0:m-1

                    FFP_view = @view FFP[comp1*nξFP+1:comp1*nξFP+nFP]

                    bc!(FFP_view, params, τ[comp1+1])

                end

                if abs(det(JFP)) > (atol + rtol * norm(xFP, Inf))

                    ΔxFP .= JFP \ FFP

                    newton_iter = 1

                    while newton_iter <= ite && norm(ΔxFP, Inf) > (atol + rtol * norm(xFP, Inf))

                        xFP .-= ΔxFP

                        update_xT!(xTFP, xFP, nξFP, nλ, nFP, m, NξFP)

                        Threads.@threads for comp1 in 1:m-1
                            xTview = @view xTFP[comp1, :]
                            taylorinteg_wrap_optim!(g!, bc!, lims, xTview, 0.0, Δτ, integtol, caches_FP[comp1], params; maxsteps=integmaxsteps)
                        end

                        FPJacobian!(JFP, xTFP, xFP, dx0FP, Δτ, nξ, nξFP, nλ, m, NξFP, NFP)
                        FPSystem!(FFP, xTFP, xFP, x0FP, dx0FP, Δτ, nξ, nξFP, m, NFP)

                        for comp1 in 0:m-1

                            FFP_view = @view FFP[comp1*(nξFP)+1:comp1*(nξFP)+nFP]

                            bc!(FFP_view, params, τ[comp1+1])

                        end

                        if abs(det(JFP)) < (atol + rtol * norm(xFP, Inf))
                            # @warn("System Jacobian was been singular. Stopping continuation.")
                            break
                        end

                        ΔxFP .= JFP \ FFP

                        newton_iter += 1

                    end

                end

                # @show xFP[1:nξ]
                # @show xFP[end-1:end]

                comp3 = 0
                for comp1 in 1:m
                    for comp2 in 1:nξ
                        comp3 += 1
                        ξFP_eval[comp1, comp2] = xFP[comp3]
                    end
                    for comp2 in 1:nξ
                        comp3 += 1
                        vFP_eval[comp1, comp2] = xFP[comp3]
                    end
                    for comp2 in 1:nλ
                        ξFP_eval[comp1, nξ+comp2] = xFP[NξFP+comp2]
                    end
                end

                # @show xFP_eval[1, :]

                EFP = H(ξFP_eval[1, :], params, 0.0)

                # @show EFP

                update_Mvec!(Mvec, xTFP, nξ, m)

                pS = pschur!(Mvec, :L)

                # @show pS.values

                push!(FP, HamiltonianFlipPoint(copy(τ), copy(ξFP_eval),
                    copy(EFP), copy(vFP_eval), Complex.(pS.values),
                    step - 1, norm(ΔxFP, Inf)))

            end

            x0 .= x1
            Φ0 .= Φ1
            E0 = E1
            FPtest0 = FPtest1
            BCPtest0 = BCPtest1

            if step % progress_step == 0 || step == maxsteps
                percentage = round(100.0 * step / maxsteps, digits=2)
                print(" \r Progress : $percentage % \t | #FP = $(length(FP)) \t |                      ")
            end

            step += 1
            Δs = Δs_ini
            astep = 0

        elseif astep < aite

            Δs = Δs * 1.e-1

            astep += 1

            if step % progress_step == 0 || step == maxsteps
                percentage = round(100.0 * step / maxsteps, digits=2)
                print(" \r Progress : $percentage % \t | #FP = $(length(FP)) \t | Adaptative step ... ")
            end

        else

            println("")
            @warn("System norm did not converge to tolerance. Stopping continuation.")
            break

        end

    end

    println(" ")

    return HamiltonianPeriodicBranch(τ,
        x[1:step-1, :, :],
        E[1:step-1],
        BCP,
        FP,
        μ[1:step-1, :],
        stbl[1:step-1],
        tol[1:step-1])

end
