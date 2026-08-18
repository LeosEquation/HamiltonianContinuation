abstract type AbstractBranchPointBranch{U<:Real} <: AbstractBranch{U} end

struct BranchPointBranch{U<:Real} <: AbstractBranchPointBranch{U}
    x::Array{U,2}                             # Solution vectors along branch
    λ::Array{Complex{U},2}               # Eigenvalues along branch
    β::Array{U,1}
    tol::Array{U,1}
end

struct HamiltonianBranchPointBranch{U<:Real} <: AbstractBranchPointBranch{U}
    x::Array{U,2}                             # Solution vectors along branch
    E::Array{U,1}                             # Energy values along branch
    λ::Array{Complex{U},2}                  # Eigenvalues along branch
    β::Array{U,1}
    tol::Array{U,1}
end

function Base.show(io::IO, obj::BranchPointBranch)
    println(io, "Branch Point Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " ├─ Parameter 1 range: [", round(minimum(obj.x[:, end-1]), digits=4), ", ", round(maximum(obj.x[:, end-1]), digits=4), "]")
    println(io, " └─ Parameter 2 range: [", round(minimum(obj.x[:, end]), digits=4), ", ", round(maximum(obj.x[:, end]), digits=4), "]")
end


function Base.show(io::IO, obj::HamiltonianBranchPointBranch)
    println(io, "Branch Point Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " ├─ Parameter 1 range: [", round(minimum(obj.x[:, end-1]), digits=4), ", ", round(maximum(obj.x[:, end-1]), digits=4), "]")
    println(io, " ├─ Parameter 2 range: [", round(minimum(obj.x[:, end]), digits=4), ", ", round(maximum(obj.x[:, end]), digits=4), "]")
    println(io, " └─ Energy range: [", round(minimum(obj.E), digits=4), ", ", round(maximum(obj.E), digits=4), "]")
end