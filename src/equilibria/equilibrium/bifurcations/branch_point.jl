struct BranchPoint{U<:Real}
    x::Vector{U}
    values::Vector{Complex{U}}
    tangent1::Vector{U}
    tangent2::Vector{U}
end