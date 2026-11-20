function Continuation(
    problem::AbstractProblem,
    config::ContinuationConfig, Δs)

    ###
    branch = build_branch(problem, config)

    ###
    test = build_test()

    ###

    continuation_point = build_continuation_point(problem, Δs)

    ###
    v = zeros(U, n)
    v[n] = one(U)

    ###
    bifurcations = build_bifurcation_points(problem)

    # bifurcations[i].type = {HamiltonianLimitPoint, HamiltonianBranchPoint}
    # bifurcations[i].point = ContinuationPoint{U}

    ###
    correction = build_correction(problem)

    # correction.regular == Array{U, n}
    # correction.bifurcation[i] == Array{U, m}

    ###
    first_step!(problem, continuation_point, branch)

    step = 2

    while step <= maxsteps

        converged = false
        newton_step = 0

        while !converged

            # Predicción del siguiente punto usando el Δs actual
            next_step!(problem, contiuation_point, Δs)

            if stop(problem, continuation_point)
                return build_branch(branch, step)
            end

            update_step!(problem, continuation_point)

            if not_regular(continuation_point, config)
                break
            end

            update_correction!(continuation_point, correction)

            newton_step = 0

            while newton_step < max_newton_ite &&
                not_convergence(continuation_point, config)

                update_step!(continuation_point, correction)

                if stop(problem, continuation_point)
                    return build_branch(branch, step)
                end

                update_step!(problem, continuation_point)

                if not_singular(continuation_point, config)
                    break
                end

                update_correction!(continuation_point, correction)

                newton_step += 1
            end

            converged = regular_convergence(continuation_point, correction, config)

            if converged
                break
            end

            # Falló Newton: reducir el paso y reintentar
            Δs *= reduce_step

            if Δs < Δs_min
                @warn("Step size became too small. Exiting.")
                return build_branch(branch, step)
            end

            previous = update_previous!(current)

        end

        if !converged
            break
        end

        # Aceptar el punto
        update_tangent!(continuation_point, v)
        update_test!(problem, continuation_point, previous)
        update_branch!(branch, continuation_point)
        # find_bifurcations!()
        update_step!(continuation_point, previous)

        step += 1

        update_step_size!(newton_step, Δs, config)

        # # Adaptar el paso según el esfuerzo de Newton
        # if newton_step <= ite_fast
        #     Δs = min(Δs * increase_step, Δs_max)
        # elseif newton_step >= ite_slow
        #     Δs = max(Δs * reduce_step, Δs_min)
        # end

    end

    return build_branch(branch, step)

end
