function EquilibriaContinuation(
    H,
    f!,
    bc!,
    g,
    x_ini::Array{U,1},
    params,
    Δs_ini::U,
    maxsteps::T;
    atol::U=eps(),
    rtol::U=eps(),
    ite::T=20,
    reduce_step::U=1.e-1,
    adapt_step::T=3,
    stbl_ord::T=1,) where {U<:Real,T<:Integer}

    ###
    n = length(x_ini)
    ordtup = [ntuple(k -> count(==(k), (i, j)), TaylorSeries.get_numvars()) for i in 1:n, j in 1:n]

    ###
    x = Array{U,2}(undef, maxsteps, n)
    E = Array{U,1}(undef, maxsteps)
    λ = Array{Complex{U},2}(undef, maxsteps, n - 1)
    ν = Array{Complex{U},3}(undef, maxsteps, n - 1, n - 1)
    lp = HamiltonianLimitPoint{U}[]
    bp = HamiltonianBranchPoint{U}[]
    stbl = Array{Int,1}(undef, maxsteps)
    tols = Array{U,1}(undef, maxsteps)

    ###
    δx = TaylorN.(1:n, order=1)
    δy = TaylorN.(1:n, order=2)
    δz = TaylorN.(1:TaylorSeries.get_numvars(), order=2)
    xaux = copy(δx)
    dx = zero(δx)
    yaux = copy(δy)
    dy = zero(δy)
    zaux = TaylorN.(1:n, order=stbl_ord)
    Haux = TaylorN(1, order=stbl_ord)

    ###
    Δs = Δs_ini
    adapt_ite = 0

    ###
    x0 = copy(x_ini)
    lptest0 = 0.0
    bptest0 = 0.0
    Φ0 = zero(x_ini)

    ###
    x1 = copy(x_ini)
    lptest1 = 0.0
    bptest1 = 0.0
    Φ1 = zero(x_ini)

    ###
    Jeval = zeros(U, n - 1, n)
    dxeval = zeros(U, n - 1)
    HessH = zeros(U, n - 1, n - 1)

    ###
    λ_ = Array{Complex{U},1}(undef, n - 1)

    ###
    y = copy(x_ini)

    ###
    J = zeros(U, n, n)
    F = zeros(U, n)
    v = zeros(U, n)
    v[n] = one(U)
    Δx = zeros(U, n)

    ###
    JLP = zeros(U, 2 * n - 1, 2 * n - 1)
    FLP = zeros(U, 2 * n - 1)
    qLP = zeros(U, 2 * n - 1)
    ΔqLP = zeros(U, 2 * n - 1)
    vLP = zeros(U, n - 1)

    ###
    JBP = zeros(U, 2 * n, 2 * n)
    FBP = zeros(U, 2 * n)
    qBP = zeros(U, 2 * n)
    ΔqBP = zeros(U, 2 * n)
    ΦBP = zeros(U, n)
    wBP = zeros(U, n - 1)
    h1 = zero(U)
    h0 = zero(U)

    ###
    x[1, :] .= x0

    for i in 1:n
        xaux[i][0][1] = x1[i]
        yaux[i][0][1] = x1[i]
        zaux[i][0][1] = x1[i]
    end

    bc!(xaux, params, 0.0)
    f!(dx, xaux, params, 0.0)

    for j in 1:n
        for i in 1:n-1
            Jeval[i, j] = dx[i][1][j]
        end
    end

    NS = nullspace(Jeval)

    if size(NS, 2) == 0
        error("Initial point does not have a valid branch direction. Cannot continue.")
    end

    if size(NS, 2) > 1
        @show NS
        error("Initial point is a codimension-2 bifurcation. Cannot continue.")
    end

    Φ0 .= NS

    EquilibriaJacobian!(J, dx, Φ0, n)
    EquilibriaSystem!(F, dx, x1, x0, Φ0, Δs, n)

    tols[1] = norm(F[1:end-1], Inf)

    for j in 1:n
        for k in 1:n-1
            Jeval[k, j] = dx[k][1][j]
        end
    end


    bc!(x1, params, 0.0)
    E[1] = H(x1, params, 0.0)

    # Haux .= H(zaux, params, 0.0)

    # for j in 1:n-1
    #     for k in 1:n-1
    #         HessH[k, j] = differentiate(ordtup[j, k], Haux[stbl_ord](δz))
    #     end
    # end

    eigeq = eigen(Jeval[:, 1:n-1])
    λ_ .= eigeq.values
    ν[1, :, :] .= eigeq.vectors
    # BiProduct!(2.0, Jeval, LinearAlgebra.I, C, n - 1)

    lptest0 = Φ0[n]
    bptest0 = det(J)

    for j in 1:n-1
        λ[1, j] = λ_[j]
    end

    # if any(real(λ_) .> (atol + rtol * norm(x1, Inf)))
    #     stbl[1] = -1
    # else
    #     if all(real(eigvals(HessH)) .< -(atol + rtol * norm(x1, Inf))) || all(real(eigvals(HessH)) .> -(atol + rtol * norm(x1, Inf)))
    #         stbl[1] = 1
    #     else
    #         stbl[1] = 0
    #     end
    # end

    if any(real(λ_) .> (atol + rtol * norm(x1, Inf)))
        stbl[1] = 0
    else
        stbl[1] = 1
    end

    #
    step = 2

    while step <= maxsteps

        for k in 1:n
            x1[k] = x0[k] + Δs * Φ0[k]
            xaux[k][0][1] = x1[k]
        end

        if g(x1, params, 0.0)
            break
        end

        bc!(xaux, params, 0.0)
        f!(dx, xaux, params, 0.0)
        EquilibriaJacobian!(J, dx, Φ0, n)
        EquilibriaSystem!(F, dx, x1, x0, Φ0, Δs, n)

        if abs(det(J)) < (atol + rtol * norm(x1, Inf))
            # @warn("System Jacobian was been singular. Stopping continuation.")
            break
        end

        Δx .= J \ F

        j = 1

        # println(" $j, $(norm(Δx, Inf)), $( atol + rtol*norm(x1, Inf) )")

        while j < ite && norm(Δx, Inf) > (atol + rtol * norm(x1, Inf))

            x1 .-= Δx

            if g(x1, params, 0.0)
                x1 .+= J \ F
                break
            end

            for k in 1:n
                xaux[k][0][1] = x1[k]
            end

            bc!(xaux, params, 0.0)
            f!(dx, xaux, params, 0.0)
            EquilibriaJacobian!(J, dx, Φ0, n)
            EquilibriaSystem!(F, dx, x1, x0, Φ0, Δs, n)

            if abs(det(J)) < (atol + rtol * norm(x1, Inf))
                # @warn("System Jacobian was been singular. Stopping continuation.")
                break
            end

            Δx .= J \ F

            j += 1

            # println(" $j, $(norm(Δx, Inf)), $( atol + rtol*norm(x1, Inf) )")

        end

        if norm(Δx, Inf) <= (atol + rtol * norm(x1, Inf))

            for k in 1:n
                zaux[k][0][1] = x1[k]
            end

            Φ1 .= J \ v

            normalize!(Φ1)

            for j in 1:n
                for k in 1:n-1
                    Jeval[k, j] = dx[k][1][j]
                end
            end

            eigeq = eigen(Jeval[:, 1:n-1])
            λ_ .= eigeq.values
            ν[step, :, :] .= eigeq.vectors
            # BiProduct!(2.0, Jeval, LinearAlgebra.I, C, n - 1)

            lptest1 = Φ1[n]
            bptest1 = det(J)

            for j in 1:n-1
                λ[step, j] = λ_[j]
            end

            # bc!(x1, params, 0.0)
            E[step] = H(x1, params, 0.0)

            # Haux .= H(zaux, params, 0.0)

            # for j in 1:n-1
            #     for k in 1:n-1
            #         HessH[k, j] = differentiate(ordtup[j, k], Haux[2](δz))
            #     end
            # end

            # if any(real(λ_) .> (atol + rtol * norm(x1, Inf)))
            #     stbl[i] = -1
            # else
            #     if all(real(eigvals(HessH)) .< -(atol + rtol * norm(x1, Inf))) || all(real(eigvals(HessH)) .> (atol + rtol * norm(x1, Inf)))
            #         stbl[i] = 1
            #     else
            #         stbl[i] = 0
            #     end
            # end

            if any(real(λ_) .> (atol + rtol * norm(x1, Inf)))
                stbl[step] = 0
            else
                stbl[step] = 1
            end

            if sign(lptest0) != sign(lptest1) &&
               sign(bptest0) == sign(bptest1) &&
               abs(lptest1 - lptest0) > (atol + rtol * norm(abs(lptest1), Inf)) &&
               step > 2

                for j in 1:n
                    y[j] = 0.5 * (x1[j] + x0[j])
                    yaux[j][0][1] = y[j]
                end

                bc!(yaux, params, 0.0)
                f!(dy, yaux, params, 0.0)

                for j in 1:n
                    for k in 1:n-1
                        Jeval[k, j] = dy[k][1][j]
                    end
                end

                Dxeigen = eigen(Jeval[:, 1:n-1])
                idx = argmin(abs.(Dxeigen.values))

                for j in 1:n-1
                    vLP[j] = real(Dxeigen.vectors[j, idx])
                end

                LPFinding!(f!, bc!, g, y, params, vLP,
                    qLP, JLP, FLP, ΔqLP,
                    dxeval, Jeval,
                    dy, yaux, ordtup, n,
                    ite, atol, rtol)

                # if norm(FLP) > biftol
                #     @warn("The limit point $(length(lp)+1) tolerance was exceded; \n norm of extended system is $(norm(FLP))")
                # end 

                # bc!(y, params, 0.0)

                push!(lp, HamiltonianLimitPoint(copy(y),
                    H(y, params, 0.0),
                    copy(Jeval),
                    copy(vLP),
                    norm(ΔqLP, Inf),
                    step - 1))

            end

            if sign(bptest0) != sign(bptest1) &&
               abs(bptest1 - bptest0) > (atol + rtol * norm(abs(bptest1), Inf)) &&
               step > 2

                for j in 1:n
                    y[j] = 0.5 * (x1[j] + x0[j])
                    yaux[j][0][1] = y[j]
                end

                bc!(yaux, params, 0.0)
                f!(dy, yaux, params, 0.0)

                for j in 1:n-1
                    for k in 1:n-1
                        Jeval[j, k] = dy[k][1][j]
                    end
                end

                Dxeigen = eigen(Jeval[:, 1:n-1])

                idx = argmin(abs.(Dxeigen.values))

                for j in 1:n-1
                    wBP[j] = real(Dxeigen.vectors[j, idx])
                end

                BPFinding!(f!, bc!, g, y, params, wBP,
                    qBP, JBP, FBP, ΔqBP,
                    dxeval, Jeval,
                    dy, yaux, ordtup, n,
                    ite, atol, rtol)

                # if norm(FBP) > biftol
                #     @warn("The branch point $(length(bp)+1) tolerance was exceded; \n norm of extended system is $(norm(FBP))")
                # end 

                h1 = sqrt(sum((x1[j] - y[j])^2 for j in 1:n))
                h0 = sqrt(sum((x0[j] - y[j])^2 for j in 1:n))

                for j in 1:n
                    ΦBP[j] = (h0 * Φ0[j] + h1 * Φ1[j]) / Δs
                end

                bc!(y, params, 0.0)

                push!(bp, HamiltonianBranchPoint(copy(y),
                    H(y, params, 0.0),
                    copy(Jeval),
                    copy(wBP),
                    copy(ΦBP),
                    norm(ΔqBP, Inf),
                    step - 1))

            end

            # bc!(x1, params, 0.0)

            x[step, :] .= x1

            tols[step] = norm(Δx, Inf)

            x0 .= x1
            Φ0 .= Φ1
            bptest0 = bptest1
            lptest0 = lptest1

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

    return HamiltonianEquilibriumBranch(x[1:step-1, :], E[1:step-1], λ[1:step-1, :], ν[1:step-1, :, :], lp, bp, stbl[1:step-1], tols[1:step-1])

end
