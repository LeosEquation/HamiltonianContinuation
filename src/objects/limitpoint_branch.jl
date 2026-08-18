abstract type AbstractLimitPointBranch{U<:Real} <: AbstractBranch{U} end

struct LimitPointBranch{U<:Real} <: AbstractLimitPointBranch{U}
    x::Array{U,2}                         # Solution vectors along branch
    λ::Array{Complex{U},2}  
    tol::Array{U, 1}                 # Eigenvalues along branch
end

struct HamiltonianLimitPointBranch{U<:Real} <: AbstractLimitPointBranch{U}
    x::Array{U,2}                             # Solution vectors along branch
    E::Array{U,1}                             # Energy values along branch
    λ::Array{Complex{U},2}
    tol::Array{U, 1}                  # Eigenvalues along branch
end

function Base.show(io::IO, obj::LimitPointBranch)
    println(io, "Limit Point Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " ├─ Parameter 1 range: [", round(minimum(obj.x[:, end-1]), digits=4), ", ", round(maximum(obj.x[:, end-1]), digits=4), "]")
    println(io, " └─ Parameter 2 range: [", round(minimum(obj.x[:, end]), digits=4), ", ", round(maximum(obj.x[:, end]), digits=4), "]")
end


function Base.show(io::IO, obj::HamiltonianLimitPointBranch)
    println(io, "Limit Point Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " ├─ Parameter 1 range: [", round(minimum(obj.x[:, end-1]), digits=4), ", ", round(maximum(obj.x[:, end-1]), digits=4), "]")
    println(io, " ├─ Parameter 2 range: [", round(minimum(obj.x[:, end]), digits=4), ", ", round(maximum(obj.x[:, end]), digits=4), "]")
    println(io, " └─ Energy range: [", round(minimum(obj.E), digits=4), ", ", round(maximum(obj.E), digits=4), "]")
end