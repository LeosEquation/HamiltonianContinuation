function update_xT!(
    xT::Array{TaylorN{U},2},
    x::Array{U,1},
    nξ::T,
    nλ::T,
    n::T,
    m::T,
    Nξ::T,) where {U<:Real,T<:Integer}

    for comp1 in eachindex(xT)
        for comp2 in 1:n
            xT[comp1][1][comp2] = 0.0
        end
    end

    comp3 = 0
    for comp1 in 1:m-1
        for comp2 in 1:nξ
            comp3 += 1
            xT[comp1, comp2][0][1] = x[comp3]
        end
        for comp2 in 1:nλ
            xT[comp1, nξ+comp2][0][1] = x[Nξ+comp2]
        end
        for comp2 in 1:n
            xT[comp1, comp2][1][comp2] = 1.0
        end
    end

end

function update_Mvec!(
    Mvec::Vector{Array{U,2}},
    xT::Array{TaylorN{U},2},
    nξ::T,
    m::T,) where {U<:Real,T<:Integer}

    for comp1 in 1:m-1
        for comp2 in 1:nξ
            for comp3 in 1:nξ
                Mvec[comp1][comp2, comp3] = xT[comp1, comp2][1][comp3]
            end
        end
    end

end