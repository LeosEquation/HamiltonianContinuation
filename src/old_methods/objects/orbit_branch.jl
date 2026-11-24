abstract type AbstractPeriodicBranch{U<:Real} <: AbstractBranch{U} end

struct PeriodicBranch{U<:Real}
    t::LinRange{U}
    x::Array{U,3}
    bcp::Vector{BranchCyclePoint{U}}
    fp::Vector{FlipPoint{U}}
    μ::Array{Complex{U},2}
    stbl::Array{Bool,1}
    tol::Array{U,1}
end

struct HamiltonianPeriodicBranch{U<:Real}
    t::LinRange{U}
    x::Array{U,3}
    E::Array{U,1}
    bcp::Vector{HamiltonianBranchCyclePoint{U}}
    fp::Vector{HamiltonianFlipPoint{U}}
    μ::Array{Complex{U},2}
    stbl::Array{Bool,1}
    tol::Array{U,1}
end

function Base.show(io::IO, obj::PeriodicBranch)
    println(io, "Periodic Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " └─ bifurcation parameter range: [", round(minimum(obj.x[:, 1, end]), digits=4), ", ", round(maximum(obj.x[:, 1, end]), digits=4), "]")

    if !isempty(obj.lcp)
        println(io, "\n  Branch Cycle Points:")
        for (i, p) in enumerate(obj.bcp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.x[1, end], digits=4), 8))
        end
    end

    if !isempty(obj.fp)
        println(io, "\n  Flip Points:")
        for (i, p) in enumerate(obj.fp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.x[1, end], digits=4), 8))
        end
    end

end

function Base.show(io::IO, obj::HamiltonianPeriodicBranch)
    println(io, "Periodic Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " └─ bifurcation parameter range: [", round(minimum(obj.E), digits=4), ", ", round(maximum(obj.E), digits=4), "]")

    if !isempty(obj.bcp)
        println(io, "\n  Branch Cycle Points:")
        for (i, p) in enumerate(obj.lcp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.E, digits=4), 8))
        end
    end

    if !isempty(obj.fp)
        println(io, "\n  Flip Points:")
        for (i, p) in enumerate(obj.fp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.E, digits=4), 8))
        end
    end

end