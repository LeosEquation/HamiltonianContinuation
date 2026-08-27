# src/integration/taylorinteg_optim.jl

import TaylorIntegration:taylorstep!

function taylorinteg_wrap_optim!(
    f!,
    bc!,
    lims,
    q0::AbstractVector{U},
    t0::T,
    tmax::T,
    abstol::T,
    cache::VectorCacheOptim,
    params;
    maxsteps::Int=500,
    reltol::T=zero(T),
) where {T<:Real,U<:Number}

    (; xaux, t, x, dx, rv, parse_eqs) = cache

    # Initial conditions
    bc!(q0, params, t0)
    if lims(q0, params, t0)
        return false
    end
    update_cache!(cache, t0, q0)
    sign_tstep = copysign(1, tmax - t0)

    # Integration
    nsteps = 1
    while sign_tstep * t0 < sign_tstep * tmax
        δt = taylorstep!(Val(parse_eqs), f!, t, x, dx, xaux, abstol, params, rv, reltol) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        if iszero(δt)
            return false
        end
        δt = sign_tstep * min(δt, sign_tstep * (tmax - t0))
        evaluate!(x, δt, q0) # new initial condition
        t0 += δt
        bc!(q0, params, t0)
        if lims(q0, params, t0)
            return false
        end
        update_cache!(cache, t0, q0)
        nsteps += 1
        if nsteps > maxsteps
            return false
        end

    end

    return true

end