struct ContinuationState{T}
    x::Vector{T}
    x0::Vector{T}
    T::Vector{T}
    T0::Vector{T}
    J::Matrix{T}
    F::Vector{T}
    Δx::Vector{T}
    v::Vector{T}
end

function build_continuation_state(C, problem)

    x0 = build_variables(C, problem)
    x = copy(x0)

    F = build_field(C, problem)
    J, T0 = build_jacobian(C, problem)

    T = copy(T0)
    Δx = J \ F

    v = zero(state.x0)
    v[end] = one(U)

    return ContinuationState(
        x,
        x0,
        T,
        T0,
        J,
        F,
        Δx
    )
end
