struct EquilibriumBifurcationPoint{U<:Real}
    dx::Vector{TaylorN{U}}
    x::Vector{TaylorN{U}}
    params::Any
    t::U
    tangent::Vector{U}
    q::Vector{U}
    F::Vector{U}
    J::Matrix{U}
end

