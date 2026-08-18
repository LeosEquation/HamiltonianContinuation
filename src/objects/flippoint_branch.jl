struct HamiltonianFlipPointBranch{U<:Real}
    t::LinRange{U}
    x::Array{U,3}
    v::Array{U,3}
    μ::Array{Complex{U},2}
    E::Array{U,1}
    tol::Array{U,1}
    step_size::Array{U,1}
end

struct HamiltonianBranchCyclePointBranch{U<:Real}
    t::LinRange{U}
    x::Array{U,3}
    v::Array{U,3}
    η::Array{U,1}
    μ::Array{Complex{U},2}
    E::Array{U,1}
    tol::Array{U,1}
    step::Array{U,1}
end

function Base.show(io::IO, obj::HamiltonianFlipPointBranch)
    println(io, "Flip Point Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " ├─ Parameter: [", minimum(obj.x[1, :, end]), ", ", maximum(obj.x[end, :, end]), "]")
    println(io, " └─ Energy : [", minimum(obj.E), ", ", maximum(obj.E), "]")
end

function Base.show(io::IO, obj::HamiltonianBranchCyclePointBranch)
    println(io, "Brancg Cycle Point Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " ├─ Parameter: [", minimum(obj.x[1, :, end]), ", ", maximum(obj.x[end, :, end]), "]")
    println(io, " └─ Energy : [", minimum(obj.E), ", ", maximum(obj.E), "]")
end