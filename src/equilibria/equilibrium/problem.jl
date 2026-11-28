struct EquilibriumProblem{U<:Real,F,B,S,P} <: AbstractProblem
    f!::F
    bc!::B
    stop::S
    params::P
    t0::U
    dim_variables::Int
    dim_params::Int
    dim::Int
    x_ini::Vector{U}
end

function EquilibriumProblem(;
    f!,
    bc!,
    stop,
    x_ini,
    params,
    t0,
    dim_params
)
    dim = length(x_ini)
    dim_variables = dim - dim_params

    return EquilibriumProblem(
        f!,
        bc!,
        stop,
        params,
        t0,
        dim_variables,
        dim_params,
        dim,
        x_ini,
    )
end