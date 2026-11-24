struct LimitPoint{U<:Real}
    x::Vector{U}
    values::Vector{Complex{U}}
    tangent::Vector{U}
end
