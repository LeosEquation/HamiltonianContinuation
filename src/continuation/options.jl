struct ContinuationOptions{U<:Real}
    atol::U
    rtol::U
    maxsteps::Int
    Δs::U
    maxite::Int
    Δsmin::U
    Δsmax::U
    target_ite::Int
    growth_factor::U
    shrink_factor::U
    max_retries::Int
end

function set_options(
    atol::U,
    rtol::U,
    maxsteps::Int,
    Δs::U,
    maxite::Int;
    Δsmin=Δs / 64,
    Δsmax=8Δs,
    target_ite=4,
    growth_factor=1.5,
    shrink_factor=0.5,
    max_retries=8
) where {U<:Real}

    return ContinuationOptions(
        atol, rtol, maxsteps, Δs, maxite,
        Δsmin, Δsmax, target_ite,
        growth_factor, shrink_factor, max_retries
    )
end