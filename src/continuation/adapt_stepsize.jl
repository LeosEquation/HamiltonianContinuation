function adapt_stepsize(Δs_current, ite, options)

    if ite <= options.target_ite - 1

        Δs_current = min(
            Δs_current * options.growth_factor,
            options.Δsmax
        )

    elseif ite >= options.target_ite + 1

        Δs_current = max(
            Δs_current * options.shrink_factor,
            options.Δsmin
        )

    end

    return Δs_current

end