abstract type AbstractBranch{U<:Real} end

abstract type AbstractEquilibriumBranch{U<:Real} <: AbstractBranch{U} end

struct EquilibriumBranch{U<:Real} <: AbstractEquilibriumBranch{U}
    x::Array{U,2}                             # Solution vectors along branch
    λ::Array{Complex{U},2}                  # Eigenvalues along branch
    lp::Vector{LimitPoint{U}}                # Limit points
    hp::Vector{HopfPoint{U}}                 # Hopf points
    bp::Vector{BranchPoint{U}}               # Branch points
    stbl::Vector{Int}                        # Stability indicators
end

struct HamiltonianEquilibriumBranch{U<:Real} <: AbstractEquilibriumBranch{U}
    x::Array{U,2}                            # Solution vectors along branch
    E::Vector{U}                             # Energy values along branch
    values::Array{Complex{U},2}              # Eigenvalues along branch
    vectors::Array{Complex{U},3}            # Eigenvectors along branch
    lp::Vector{HamiltonianLimitPoint{U}}     # Limit points
    bp::Vector{HamiltonianBranchPoint{U}}    # Branch points
    stability::Vector{Int}                        # Stability indicators
    tol::Array{U,1}                         # Tolerance in each step
end

function Base.show(io::IO, obj::EquilibriumBranch)
    println(io, "Equilibrium Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " └─ bifurcation parameter range: [", round(minimum(obj.x[:, end]), digits=4), ", ", round(maximum(obj.x[:, end]), digits=4), "]")

    if !isempty(obj.lp)
        println(io, "\n  Limit Points:")
        for (i, p) in enumerate(obj.lp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.x[end], digits=4), 8))
        end
    end

    if !isempty(obj.bp)
        println(io, "\n  Branch Points:")
        for (i, p) in enumerate(obj.bp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.x[end], digits=4), 8))
        end
    end

    if !isempty(obj.hp)
        println(io, "\n  Hopf Points:")
        for (i, p) in enumerate(obj.hp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.x[end], digits=4), 8))
        end
    end
end

function Base.show(io::IO, obj::HamiltonianEquilibriumBranch)
    println(io, "Hamiltonian Equilibrium Branch")
    println(io, " ├─ # steps:        ", size(obj.x, 1))
    println(io, " ├─ parameter range: [", round(minimum(obj.x[:, end]), digits=4), ", ", round(maximum(obj.x[:, end]), digits=4), "]")
    println(io, " └─ energy range:   [", round(minimum(obj.E), digits=4), ", ", round(maximum(obj.E), digits=4), "]")

    if !isempty(obj.lp)
        println(io, "\n  Hamiltonian Limit Points:")
        for (i, p) in enumerate(obj.lp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.x[end], digits=4), 8), " |  E: ", round(p.E, digits=4))
        end
    end

    if !isempty(obj.bp)
        println(io, "\n  Hamiltonian Branch Points:")
        for (i, p) in enumerate(obj.bp)
            println(io, "   (", i, ") idx: ", lpad(p.idx, 4), "  |  bif param: ", rpad(round(p.x[end], digits=4), 8), " |  E: ", round(p.E, digits=4))
        end
    end
end

struct HamiltonianEquilibriumPoint{U<:Real}
    x::Vector{U}
    E::U
    values::Vector{Complex{U}}
    vectors::Matrix{Complex{U}}
    stability::Int
    tol::U
end

function Base.show(io::IO, ::MIME"text/plain", p::HamiltonianEquilibriumPoint{U}) where {U}

    println(io, "Hamiltonian Equilibrium Point")
    println(io, "────────────────────────────")

    println(io, "Parameter : ", p.x[end])
    println(io, "Stability : ", p.stability)
    @printf(io, "Tolerance : %.2e\n", p.tol)
    println(io)

    println(io, "Eigenpairs:")
    println(io)

    for i in eachindex(p.values)

        λ = p.values[i]
        v = p.vectors[:, i]

        println(io, "  λ$i = $λ")
        println(io, "     v$i = [")

        for comp in v
            println(io, "        $comp")
        end

        println(io, "     ]")
        println(io)
    end
end

function Base.getindex(eqb::HamiltonianEquilibriumBranch{U}, i::Int) where {U}
    HamiltonianEquilibriumPoint{U}(
        eqb.x[i, :],
        eqb.E[i],
        eqb.values[i, :],
        eqb.vectors[i, :, :],
        eqb.stability[i],
        eqb.tol[i]
    )
end

Base.length(eqb::HamiltonianEquilibriumBranch) = length(eqb.E)

Base.size(eqb::HamiltonianEquilibriumBranch) = (length(eqb),)

