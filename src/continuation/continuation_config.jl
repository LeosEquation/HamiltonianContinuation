struct ContinuationConfig{T<:Integer,U<:Real}
    maxsteps::T
    atol::U
    rtol::U
    max_newton_ite::T
    ite_fast::T
    ite_slow::T
    increase_step::U
    reduce_step::U
    Δs_min::U
    Δs_max::U
end

function build_config(maxsteps, atol, rtol, Δs)
    return ContinuationConfig(
        maxsteps,
        atol,
        rtol,
        20,              # max_newton_ite
        4,               # ite_fast
        10,              # ite_slow
        1.25,            # increase_step
        0.5,             # reduce_step
        0.01 * abs(Δs),  # Δs_min
        5.0 * abs(Δs)    # Δs_max
    )
end