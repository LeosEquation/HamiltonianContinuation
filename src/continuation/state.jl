"""
    ContinuationState{T}

Stores the numerical variables and auxiliary arrays required during
a continuation process.

# Fields
- `x::Vector{T}`: Current continuation point.
- `x0::Vector{T}`: Previous accepted continuation point.
- `T::Vector{T}`: Current tangent vector.
- `T0::Vector{T}`: Tangent vector at the previous accepted point.
- `J::Matrix{T}`: Jacobian matrix of the system.
- `F::Vector{T}`: Vector field evaluated at the current point.
- `J_factor::Matrix{T}`: Auxiliary matrix used to store the LU factorization
  of the Jacobian without modifying `J`.
- `Δx::Vector{T}`: Newton correction vector.
- `v::Vector{T}`: Auxiliary vector used by the continuation algorithm.
"""
struct ContinuationState{T}
    x::Vector{T}
    x0::Vector{T}
    T::Vector{T}
    T0::Vector{T}
    J::Matrix{T}
    F::Vector{T}
    J_factor::Matrix{T}
    Δx::Vector{T}
    v::Vector{T}
end

"""
    build_continuation_state(C, problem)

Builds the numerical state required by the continuation method `C`.

Initializes the continuation variables, vector field, Jacobian, tangent
vector, and auxiliary arrays using the corresponding construction
functions.

The auxiliary matrix `J_factor` is initialized as a copy of the Jacobian,
while the Newton correction vector and auxiliary vector are allocated
without initializing their entries.

# Arguments
- `C::Type{<:AbstractContinuation}`: Continuation method to be used.
- `problem::AbstractProblem`: Problem defining the system to be continued.

# Returns
A `ContinuationState` containing the initialized numerical variables
and auxiliary arrays.
"""
function build_continuation_state(C, taylor, problem)

    x0 = build_variables(C, problem)
    x = copy(x0)

    F = build_field(C, taylor, problem)
    J, T0 = build_jacobian(C, taylor, problem)

    J_factor = copy(J)

    T = copy(T0)
    Δx = similar(F)

    v = zero(x0)
    v[end] = one(eltype(v))

    return ContinuationState(
        x,
        x0,
        T,
        T0,
        J,
        F,
        J_factor,
        Δx,
        v
    )

end