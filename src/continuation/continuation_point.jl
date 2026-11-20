struct EquilibriumContinuationPoint{U<:Real}
    dx::Vector{TaylorN{U}}
    x::Vector{TaylorN{U}}
    params::Any
    t::U
    tangent::Vector{U}
    tangent::Vector{U}
    F::Vector{U}
    J::Vector{U}
    x0::Vector{U}
    dx0::Vector{U}
    tangent0::Vector{U}
end

struct EquilibriumBifurcationPoint{U<:Real}
    dx::Vector{TaylorN{U}}
    x::Vector{TaylorN{U}}
    params::Any
    t::U
    tangent::Vector{U}
    q::Vector{U}
    F::Vector{U}
    J::Vector{U}
end





